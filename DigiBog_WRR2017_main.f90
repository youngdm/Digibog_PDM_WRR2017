PROGRAM DIGIBOG_BB

! Version 2015_03_17 - this version combines DigiBog_16_01_2013.f90
! and DigiBog_Hydro_NRW_main_v2.f90 (version 28/08/13) to model 2.5D space as
! a blanket bog to test management impacts.

!    Copyright (C) 2015 A.J. Baird, P.J. Morris, L.R. Belyea, and D.M. Young

!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.

! -----------------------------------------------------------------------------
! Section 1.0 Program header
! -----------------------------------------------------------------------------

  !Description
  !A model to simulate water tables  and peat accumulation / wastage in
  !ombrotrophic bogs and other shallow aquifers using the Boussinesq equation.
  !The model uses centimetres and seconds.
  !Details of how the model works may be obtained from the code owners
  !(email addresses below).

  !Current code owners
  !Paul J. Morris, Andy J. Baird*, Lisa R. Belyea and Dylan M. Young
  !*School of Geography
  !University of Leeds
  !Leeds
  !LS2 9JT
  !a.j.baird@leeds.ac.uk
  !l.belyea@qmul.ac.uk
  !p.j.morris@reading.ac.uk
  !d.m.young@leeds.ac.uk

  !Modification history of code
  !Programmer           Date           Modifications
  !============================================================================
  !Andrew J. Baird      08/04/2005     Original 1.5-d code
  !----------------------------------------------------------------------------
  !Paul J. Morris       24/02/2006     Conversion to 2.5-d and Fortran
  !                                    Standards implemented
  !----------------------------------------------------------------------------
  !Paul J. Morris       27/02/2006     Testing, minor corrections
  !----------------------------------------------------------------------------
  !Paul J. Morris       05/03/2006     Subroutines 'column_activation' and
  !                                    'steady_state_check' written. Former
  !                                    no longer exists (removed 18/07/2013).
  !----------------------------------------------------------------------------
  !Paul J. Morris       19/06/2006     Testing completed
  !----------------------------------------------------------------------------
  !Paul J. Morris       20/03/2007     Code cleaning
  !----------------------------------------------------------------------------
  !Paul J. Morris       24/04/2007     Further cleaning, including removal of
  !                                    replicated spatial step reference in
  !                                    'move_water' subroutine
  !----------------------------------------------------------------------------
  !Paul J. Morris       09/05/2007     Above-ground storage facilitated in
  !                                    2.5-d version
  !----------------------------------------------------------------------------
  !Paul J. Morris       02/07/2008     Final code cleaning, full annotation
  !----------------------------------------------------------------------------
  !Andy J. Baird        18/07/2013     Addition of variable net rainfall read
  !                                    in from file. Change to how boundary
  !                                    conditions specified; zero-flow
  !                                    Neumann condition also now allowed for,
  !                                    as are internal boundaries. Change to
  !                                    how output written to file (can now
  !                                    write to file multiple times before end
  !                                    of run).
  !----------------------------------------------------------------------------
  !Andy J. Baird        27/08/2013     Code cleaning and de-bugging of version
  !                                    from 18/07/2013
  !----------------------------------------------------------------------------
  !Dylan M. Young       20/05/14       Combine 1D and Hydro models for use in
  !                                    blanket bog simulation of management
  !                                    impacts. For testing, set fixed Dirichlet
  !                                    condition of column elevation * 0.5.
  !                                    Added memory allocation for x_flux and
  !                                    y_flux in wat_k_mean (Hydro_procs) for
  !                                    de-bugging purposes.
  !-----------------------------------------------------------------------------
  !Dylan M. Young       11/06/14       Reorganisation of code as suggested by
  !                                    AJB and agreed with AJB and PJM.
  !                                    Removal of de-bugging code and resetting
  !                                    x_flux and y_flux to local arrays.
  !                                    Removal of steady-state check.
  !-----------------------------------------------------------------------------
  !Dylan M. Young       18/06/14       Addition of variable Dirichlet condition
  !                                    water-table read in from the parameter
  !                                    file.
  !-----------------------------------------------------------------------------
  !Dylan M. Young       13/07/14       Added timestep multiplier for calculated
  !                                    timestep for stability improvements
  !-----------------------------------------------------------------------------
  !Dylan M. Young       01/01/15       Added mineral layer for sloping model 
  !                                    stability. Can be used for both raised 
  !                                    and blanket bog models.
  !-----------------------------------------------------------------------------
  !Dylan M. Young       17/03/15       Weather now read-in on a monthly basis
  !                                    (can be used for stable or variable 
  !                                    weather). Mean annual temp. is used for 
  !                                    new litter production.
  !-----------------------------------------------------------------------------
  !Dylan M. Young       20/03/15       Recalcitrance added to anoxic
  !                                    decomposition.
  !-----------------------------------------------------------------------------
  !Dylan M. Young       28/04/15       Added write-out for K profile and maximum 
  !                                    timestep calc.
  !-----------------------------------------------------------------------------
  !Code description.
  !Language: Fortran 95.
  !Software standards: Andrews, P. (1998). Unified Model Documentation Paper No.
  !3, Software Standards for the Unified Model: Fortran and Unix, Version 7.2,
  !Met Office, UK, 17 pp.
  !-----------------------------------------------------------------------------

  !Modules used:
  USE hydro_procedures_bb_transmax
  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Section 2.0 Definitions; Declarations
!-------------------------------------------------------------------------------

! Single-valued variable (and model parameter) declarations.

  INTEGER :: alloc_error,        &  !Memory allocation error flag
             t_extent,           &  !Number of simulated model years
             output_interval,    &  !Interval (years) for writing to output
                                    !files
             sub_year_counter,   &  !Counter for use in main model loop
             sub_year_counter_su,&  !Counter for use in main model loop
             output_counter,     &  !Output counter for use in main time loop
             year_counter,       &  !Counter for total number of model years
             week_counter,      &  !Counter for total months
             new_week,          &  !Flag for timing months
             no_layers,          &  !Number of layers
             x_extent,           &  !Model grid extent in x direction (number)
             y_extent,           &  !Model grid extent in y direction (number)
             z_extent,           &  !Spatial index counter
             x,                  &  !Spatial index counter
             y,                  &  !Spatial index counter
             z                      !Spatial index counter

  REAL(kind=q) :: oxic_decay_base,     &  !Proportion (yr^-1)
                  anoxic_decay_base,   &  !Proportion (yr^-1)
                  base_temp,           &  !deg. Celsius
                  Q10_oxic,            &  !Factor (dimensionless)
                  Q10_anoxic,          &  !Factor (dimensionless)
                  oxic_decay,          &  !Proportion (yr^-1)
                  anoxic_decay,        &  !Proportion (yr^-1)
                  density,             &  !Peat bulk density (g cm^-3)
                  porosity,            &  !Drainable porosity (cm^3 cm^-3 --
                                          !a proportion)
                  k_param_a,           &  !The a parameter in the k equation
                                          !(cm yr^-1)
                  k_param_b,           &  !The b parameter in the k equation
                                          !(dimensionless)
                  timestep,            &  !yr
                  transmax,            &  !Maximum transmissivity for timestep
                  t_step_multi,        &  !timestep multiplication factor
                  t_step_sum,          &  !sum of timesteps
                  t_step_sum_su,       &  !sum of summer timesteps
                  mean_t_step,         &  !Average annual timesteps
                  annual_tsteps,       &  !Number of annual timesteps
                  week_tsteps,        &  !Number of monthly timesteps
                  max_tstep,           &  !Max timestep value 
                  rainfall,            &  !Net rainfall (cm yr^-1)
                  temperature,         &  !deg. C
                  spatial_step,        &  !Model horizontal increment
                                          !(x & y) (cm)
                  pond_depth,          &  !Depth of surface ponding (cm)
                  mean_temp               !Mean annual temperature deg. C

 !Input and output files
 CHARACTER(LEN=38) :: data_file_name_010, &  !IN Model run information file
                      data_file_name_020, &  !IN net rainfall values
                      data_file_name_030, &  !IN temperature values
                      data_file_name_040, &  !IN column activation status file
                      data_file_name_050, &  !IN base altitude input file
                      data_file_name_060, &  !OUT Layer mass remaining output
                      data_file_name_070, &  !OUT Column height output file
                      data_file_name_080, &  !OUT Water-table height output file
                      data_file_name_090, &  !OUT Transmissivity output file
                      data_file_name_100, &  !OUT Mass per area output file
                      data_file_name_110, &  !OUT Mean annual column water table
                      data_file_name_120, &  !OUT Layer wet proportion file
                      data_file_name_130, &  !OUT timestep
                      data_file_name_140, &  !OUT timestep sum
                      data_file_name_150, &  !OUT K_profile
                      data_file_name_160     !OUT Summer water table

  REAL(KIND=q), ALLOCATABLE, DIMENSION(:,:) :: base_altitude,     & !Above datum
                                               water_change,      &
                                               water_table,       & !Above base
                                               wk_mean,           & !Depth-av. K
                                               col_wt_depth,      & !cm
                                               col_wt_sum,        & !cm
                                               col_wt_sum_su,     & !cm
                                               col_wt_depth_av,   & !cm
                                               col_wt_depth_av_su,& !cm
                                               col_mass_per_area    !g cm^-2

  REAL(KIND=q), ALLOCATABLE, DIMENSION(:,:,:) :: wet_proportion !proportion

  !layer_attributes stores layer thickness, k (cm yr^-1) and s (cm^3 cm^-3)
  !transmissivity stores layer elevation above base and transmissivity
  !layer_storage stores layer elevation (cm) above base and each layer's water
  !capacity as volume per unit area; i.e. expressed as a depth (cm)
  !layer_mass stores layer current (g cm^-2), initial (g cm^-2)
  !and remaining mass (proportion).
  REAL(KIND=q), ALLOCATABLE, DIMENSION(:,:,:,:) :: layer_attributes, &
                                                   transmissivity, &
                                                   layer_storage, &
                                                   layer_mass

  !activation_status indicates whether/how column participates in simulation
  CHARACTER(8), ALLOCATABLE, DIMENSION(:,:) :: activation_status

  CHARACTER(LEN=1) :: accept_licence !Flag to accept software licence

!-------------------------------------------------------------------------------
!Licence agreement
!-------------------------------------------------------------------------------
  WRITE (*, *)
  WRITE (*, *) "---------------------------------------------------------------"
  WRITE (*, *) "DigiBog. Copyright (C) 2015 A.J. Baird, L.R. Belyea,  &
          P.J. Morris, and D.M. Young"
  WRITE (*, *) " "
  WRITE (*, *) "This program comes with ABSOLUTELY NO WARRANTY;"
  WRITE (*, *) "This is free software, and you are welcome to redistribute it"
  WRITE (*, *) "under certain conditions; refer to accompanying licence or see"
  WRITE (*, *) "<http://www.gnu.org/licenses/>"
  WRITE (*, *)
  WRITE (*, *)
!-------------------------------------------------------------------------------
! Section 3.0 Data Input; Memory Management
!-------------------------------------------------------------------------------

  ! Name input and output data files.

  !Inputs
  data_file_name_010 = "010_DigiBog_BB_IN_information.txt"
  data_file_name_020 = "020_DigiBog_BB_IN_net_rain.txt"
  data_file_name_030 = "030_DigiBog_BB_IN_temp.txt"
  data_file_name_040 = "040_DigiBog_BB_IN_column_status.txt"
  data_file_name_050 = "050_DigiBog_BB_IN_baltitude.txt"

  !Outputs
  data_file_name_060 = "060_DigiBog_BB_OUT_layer_mass.txt"
  data_file_name_070 = "070_DigiBog_BB_OUT_column_height.txt"
  data_file_name_080 = "080_DigiBog_BB_OUT_wt_height.txt"
  data_file_name_090 = "090_DigiBog_BB_OUT_transmiss.txt"
  data_file_name_100 = "100_DigiBog_BB_OUT_mass_area.txt"
  data_file_name_110 = "110_DigiBog_BB_OUT_wt_depth.txt"
  data_file_name_120 = "120_DigiBog_BB_OUT_layer_wet_prop.txt"
  data_file_name_130 = "130_DigiBog_BB_OUT_col_t_step.txt"
  data_file_name_140 = "140_DigiBog_BB_OUT_t_step_sum.txt"
  data_file_name_150 = "150_DigiBog_BB_OUT_k_profile.txt"
  data_file_name_160 = "160_DigiBog_BB_OUT_wt_depth_summer.txt"

! Open files for input and output of data.
  OPEN(UNIT=010, FILE=data_file_name_010, STATUS="OLD")
  OPEN(UNIT=020, FILE=data_file_name_020, STATUS="OLD")
  OPEN(UNIT=030, FILE=data_file_name_030, STATUS="OLD")
  OPEN(UNIT=040, FILE=data_file_name_040, STATUS="OLD")
  OPEN(UNIT=050, FILE=data_file_name_050, STATUS="OLD")
  OPEN(UNIT=060, FILE=data_file_name_060, STATUS="REPLACE")
  OPEN(UNIT=070, FILE=data_file_name_070, STATUS="REPLACE")
  OPEN(UNIT=080, FILE=data_file_name_080, STATUS="REPLACE")
  OPEN(UNIT=090, FILE=data_file_name_090, STATUS="REPLACE")
  OPEN(UNIT=100, FILE=data_file_name_100, STATUS="REPLACE")
  OPEN(UNIT=110, FILE=data_file_name_110, STATUS="REPLACE")
  OPEN(UNIT=120, FILE=data_file_name_120, STATUS="REPLACE")
  OPEN(UNIT=130, FILE=data_file_name_130, STATUS="REPLACE")
  OPEN(UNIT=140, FILE=data_file_name_140, STATUS="REPLACE")
  OPEN(UNIT=150, FILE=data_file_name_150, STATUS="REPLACE")
  OPEN(UNIT=160, FILE=data_file_name_160, STATUS="REPLACE")

! Read data from information input data file.
  READ (010, *) oxic_decay_base
  READ (010, *) anoxic_decay_base
  READ (010, *) base_temp
  READ (010, *) Q10_oxic
  READ (010, *) Q10_anoxic
  READ (010, *) density
  READ (010, *) porosity
  READ (010, *) k_param_a
  READ (010, *) k_param_b
  READ (010, *) t_extent
  READ (010, *) annual_tsteps
  READ (010, *) output_interval
  READ (010, *) x_extent
  READ (010, *) y_extent
  READ (010, *) spatial_step
  READ (010, *) pond_depth
  READ (010, *) mean_temp
  READ (010, *) t_step_multi
  READ (010, *) max_tstep

  !Set the extent of layers
  z_extent = t_extent + 2 !  Ponding layer + added mineral + peat layers

  !Allocate memory to arrays.

  !Columns

  ALLOCATE(col_wt_depth(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*,*) "Model could not allocate space for col_wt_depth"
    STOP
  END IF

  ALLOCATE(col_wt_depth_av(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*,*) "Model could not allocate space for col_wt_depth_av"
    STOP
  END IF

  ALLOCATE(col_wt_depth_av_su(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*,*) "Model could not allocate space for col_wt_depth_av"
    STOP
  END IF

  ALLOCATE(col_wt_sum(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*,*) "Model could not allocate space for col_wt_sum"
    STOP
  END IF

  ALLOCATE(col_wt_sum_su(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*,*) "Model could not allocate space for col_wt_sum_su"
    STOP
  END IF

  ALLOCATE(col_mass_per_area(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*,*) "Model could not allocate space for col_mass_per_area"
    STOP
  END IF

  ALLOCATE(base_altitude(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for base_altitude"
    STOP
  END IF

  ALLOCATE(water_change(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for water_change"
    STOP
  END IF

  ALLOCATE(water_table(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for water_table"
    STOP
  END IF

  ALLOCATE(wk_mean(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for wk_mean"
    STOP
  END IF

  ALLOCATE(activation_status(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for activation_status"
    STOP
  END IF

  ALLOCATE(wet_proportion(x_extent, y_extent, z_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*,*) "Model could not allocate space for wet_proportion"
    STOP
  END IF

  ALLOCATE(layer_mass(x_extent, y_extent, z_extent, 3), STAT = alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for layer_mass"
    STOP
  END IF

  ALLOCATE(layer_attributes(x_extent, y_extent, z_extent, 3), &
                                                             STAT = alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for layer_attributes"
    STOP
  END IF

  ALLOCATE(transmissivity(x_extent, y_extent, z_extent, 2), STAT = alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for transmissivity"
    STOP
  END IF

  ALLOCATE(layer_storage(x_extent, y_extent, z_extent, 2), STAT = alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for layer_storage"
    STOP
  END IF


!-------------------------------------------------------------------------------
! Section 4.0 Initialisation
!-------------------------------------------------------------------------------

  !Initialise all single-valued variables and arrays.

  !Single-valued variables.
  year_counter = 1
  week_counter = 0
  annual_tsteps = 1.0
  week_tsteps = annual_tsteps / REAL(52)
  timestep = week_tsteps
  t_step_sum = 0.0
  t_step_sum_su = 0.0
  max_tstep = REAL(max_tstep) / REAL(525600)
  mean_t_step = 0.0
  sub_year_counter = 0
  sub_year_counter_su = 0
  output_counter = 0
  oxic_decay = 0.0
  anoxic_decay = 0.0
  x = 1
  y = 1
  z = 1
  no_layers = 3
  new_week = 1

  !Arrays.
  layer_mass = 0.0
  wet_proportion = 0.0
  col_wt_depth = 0.0
  col_wt_sum = 0.0
  col_wt_sum_su = 0.0
  col_wt_depth_av = 0.0
  col_wt_depth_av_su = 0.0
  col_mass_per_area = 0.0
  layer_attributes = 0.0
  water_change = 1.0
  wk_mean = 1.0
  transmissivity = 1.0
  layer_storage = 1.0
  water_table = 0.0
  base_altitude = 0.0

  !Read data from files to rank-two arrays for active columns
  DO x = 1, x_extent
    DO y = 1, y_extent
      READ (040, *) activation_status(x, y)
      READ (050, *) base_altitude(x, y)
    END DO
  END DO

  !Initialise pond properties. k and s only set for layer 1 for
  !"on" and "diri" columns because additional layers added each year.
  DO x = 1, x_extent
    DO y = 1, y_extent
      IF (activation_status(x, y) == "on") THEN
      !Pond properties
      layer_attributes(x, y, no_layers, 1) = pond_depth !thickness
      layer_attributes(x, y, no_layers, 2) = 0.0 !k
      layer_attributes(x, y, no_layers, 3) = 1.0 !s
      END IF
    END DO
  END DO

  ! Initialise the properties for a layer of mineral soil.
  DO x = 1, x_extent
    DO y = 1, y_extent
      IF (activation_status(x, y) == "on") THEN
        !Specify layer thickness of 2cm
        layer_attributes(x, y, 1, 1) = 20.0
        !Calculate layer 1 k
        layer_attributes(x, y, 1, 2) = 3153.6 !Layer K = 10-4cm s-1
        !Layer porosity
        layer_attributes(x, y, 1, 3) =  0.3 
        !soils
        !Column elevation for layer 1 = layer thickness of layer 1
        transmissivity(x, y, 1, 1) = layer_attributes(x, y, 1, 1)
        !Calculate layer 1 transmissivity = layer 1 thickness * layer 1 k
        transmissivity(x, y, 1, 2) =  transmissivity(x, y, 1, 1) &
                                       * layer_attributes(x, y, 1, 2)
        !Set wet proportion
        wet_proportion(x, y, 1) = 1.0
        !Initial column properties
        !Assign mass per area = layer 1 mass
        !col_mass_per_area(x, y) = layer_mass(x, y, 1, 1)
        !Set the water table to the top of layer 1 for active columns
        water_table(x, y) = transmissivity(x, y, 1, 1)
      END IF
    END DO
  END DO

  ! Initialise properties of first layer of peat profile for "on" columns.
  DO x = 1, x_extent
    DO y = 1, y_extent
      IF (activation_status(x, y) == "on") THEN
        !Calculations and assignments for columns with either active status or
        !Calculate layer 2 mass
        ! Equation based on Belyea and Clymo (2001) - see main model loop.
        ! Assumed that water table at surface and temp. at base value
        layer_mass(x ,y, 2, 1)  = (9.3**2.0) * 0.0001
        !layer mass = 0.008649 g cm^2
        !Assign layer 2 initial mass to be layer 2 mass
        layer_mass(x, y, 2, 2) = layer_mass(x, y, 2, 1)
        !Set layer 2 remaining mass to be 100%
        layer_mass(x, y, 2, 3) = 1.0
        !Calculate layer thickness
        layer_attributes(x, y, 2, 1) = layer_mass(x, y, 2, 1) / density
        !Calculate layer 2 k
        !Equation below from Morris et al. (2011) - see main model loop.
        layer_attributes(x, y, 2, 2) =  k_param_a * (exp (k_param_b))
        !Layer porosity
        layer_attributes(x, y, 2, 3) =  porosity
        !Column elevation for layer 2 = layer thickness of layer 1 + layer 2
        transmissivity(x, y, 2, 1) = transmissivity(x, y, 1, 1) + &
                                       layer_attributes(x, y, 2, 1)
        !Calculate layer 2 transmissivity = layer 1 thickness * layer 1 k
        transmissivity(x, y, 2, 2) =  transmissivity(x, y, 1, 1) &
                                      + (layer_attributes(x, y , 2 ,1) &
                                       * layer_attributes(x, y, 2, 2))
        !Set wet proportion
        wet_proportion(x, y, 2) = 1.0
        !Initial column properties
        !Assign mass per area = layer 2 mass (does not decompose z1 mineral)
        col_mass_per_area(x, y) = layer_mass(x, y, 2, 1)
        !Set the water table to the top of layer 2 for active columns
        water_table(x, y) = transmissivity(x, y, 2, 1)
      END IF
    END DO
  END DO

  !Set water-table height for Dirichlet condition. x1 and x_extent are
  !"neu" or "off" cells as is y_extent
  Diri_wt: DO x = 2, (x_extent -1)
    DO y = 1, (y_extent -1)
      IF (activation_status(x, y) == "diri") THEN
         water_table(x, y) =  transmissivity(x, y + 1, 2, 1) 
      END IF
    END DO
  END DO Diri_wt

!-------------------------------------------------------------------------------
! Section 5.0 Main Calculations; Time Management
!-------------------------------------------------------------------------------

  !1.Commence main annual loop.
  Annual_loop: DO

    !Display model clock on screen in units of years.
    WRITE (*,*) "Model year ", year_counter

    !2.Start sub annual timesteps
    Sub_annual_loop: DO


  !2.1 Weather -----------------------------------------------------------------
  !Check if the sum of timesteps exceed the number of monthly timesteps. 
  !Or if a new year has been signalled with new_year.
  !If so read in the next values for temperature and net rainfall
  IF (new_week == 1) THEN
    !Read next net rainfall and temperature values
    READ (020, *) rainfall
    READ (030, *) temperature
    !Calculate decay rate variables for current month.
    oxic_decay = oxic_decay_base & 
      * Q10_oxic**((temperature - base_temp) / 10.0)
    anoxic_decay = anoxic_decay_base &
      * Q10_anoxic**((temperature - base_temp) / 10.0)
    !Reset time counters for next trip through monthly loop
    WRITE(140, '(20F20.8)') t_step_sum
    t_step_sum = 0.0
    t_step_sum_su = 0.0
    new_week = 0
  ENDIF
  !----------------------------------------------------------------------------

  !Count the number of sub-annual loops (≡ annual_tsteps)
  !the count is used in the average water-table calc in the annual loop
  sub_year_counter = sub_year_counter + 1
  IF (week_counter >= 17 .AND. week_counter <= 31) THEN
    sub_year_counter_su = sub_year_counter_su + 1
  END IF

  !Re-initialise column properties
  col_mass_per_area = 0.0

  !Set the timestep for this iteration to be equal to 1 month and reset
  !transmax to 0.0. 
  !minute
  timestep = week_tsteps
  transmax = 0.0


      !Calculate amount of water (expressed as a depth) which may be stored in
      !each layer within each column
      CALL layer_water_depth(x_extent, y_extent, no_layers, layer_attributes, &
                          layer_storage, activation_status)

      !Calculate the depth-averaged K below the water table for each column
      CALL wat_k_mean(x_extent, y_extent, no_layers, transmax, water_table, &
                      wk_mean, layer_attributes, transmissivity, &
                      activation_status)

      !Reduce the calculated timestep for use with both hydrological and
      !decomposition routines
      timestep = (0.5 * ((0.5 * spatial_step) ** 2) * porosity / transmax) &
        * 0.8

      timestep = timestep * t_step_multi

      !Reduce the calculated timestep for use with both hydrological and
      !decomposition routines. Timestep set to 10 minutes for first 50 years
      !and then to a maxumum of 50 minutes thereafter
      IF (year_counter > 50) THEN
        IF (timestep > REAL(50) / REAL(525600)) THEN
          timestep = REAL(50) / REAL(525600)
        END IF
        ELSE IF(year_counter <= 50) THEN
          IF (timestep >  max_tstep) THEN
            timestep =  max_tstep
          END IF
      END IF      
      
  !Keep the sum of sub-weekly timesteps equal to the weekly timestep to
  !prevent the over-accumulation of time in the model time extent
  IF (timestep + t_step_sum > week_tsteps) THEN
    timestep = week_tsteps - t_step_sum
    !Set the new month flag
    new_week = 1
    !Add a month to the month counter
    week_counter = week_counter + 1
  END IF

  !Sum up the value of each timestep to keep track of sub-monthly looping
  t_step_sum = t_step_sum + timestep
  IF (week_counter >= 17 .AND. week_counter <= 31) THEN
    t_step_sum_su = t_step_sum_su + timestep
  END IF
  !----------------------------------------------------------------------------

      !Calculate the amount of water (expressed as a depth) that moves between
      !each column. The flow law is based on the Boussinesq equation.
      CALL move_water(x_extent, y_extent, timestep, spatial_step, rainfall, &
                      base_altitude, water_change, water_table, wk_mean, &
                      activation_status)

      !Update position of the water table in each column
      CALL water_table_update(x_extent, y_extent, no_layers, water_change, &
                             water_table, layer_attributes, layer_storage, &
                             activation_status)

      !Calculate new water-table depth for each active column and use to
      !calculate mean water-table depth for year. The mean value is used to
      !calculate the new annual layer litter production.
      Water_table_depth: DO x = 1, x_extent
        DO y = 1, y_extent
          IF (activation_status(x, y) == "on") THEN
            col_wt_depth(x, y) = transmissivity(x, y, (no_layers - 1), 1)&
                                  - water_table(x, y)
            col_wt_sum(x, y)  = col_wt_sum(x, y) + col_wt_depth(x, y)
            IF (week_counter >= 17 .AND. week_counter <= 31) THEN
                col_wt_sum_su(x, y)  = col_wt_sum_su(x, y) + col_wt_depth(x, y)
            END IF
          END IF
        END DO
      END DO Water_table_depth

      !Decompose peat in each layer using info. from previous timestep.
      !Use proportionate decay model. Amount of peat lost depends on whether
      !peat above or below the water table. Scan upwards through all model
      !layers. Use new layer properties to update column properties.
      Decompose: DO x = 1, x_extent
        DO y = 1, y_extent
          IF (activation_status(x, y) == "on") THEN
            !Decompose peat layers only
            DO z = 2, (no_layers - 1)
            !Update layer mass including a recalcitrance function (see PM)
            !applied only to anoxic decomposition
            layer_mass(x, y, z, 1) = (layer_mass(x, y, z, 1) *&
                                      wet_proportion(x, y, z) *&
                                      (exp (- anoxic_decay * timestep *&
                                      layer_mass(x, y, z, 3)))) + &
                                      (layer_mass(x, y, z, 1) *&
                                      (1 - wet_proportion(x, y, z)) *&
                                      (exp (- oxic_decay * timestep))) 
            !Update layer remaining mass
            layer_mass(x, y, z, 3) = layer_mass(x, y, z, 1) /&
                                     layer_mass(x, y, z, 2)
            !Layer thickness
            layer_attributes(x, y, z, 1) = layer_mass(x, y, z, 1) / density
            !Current layer mass for addition to mass per area array
            col_mass_per_area(x, y) = col_mass_per_area(x, y) +&
              layer_mass(x, y, z, 1)
            END DO
          END IF
        END DO
      END DO Decompose ! End of decomposition loop

      !Recalculate layer k
      !For peat layers only
      Layer_k: DO x = 1, x_extent
        DO y = 1, y_extent
          IF (activation_status(x, y) == "on") THEN
            DO z = 2, (no_layers - 1)
            layer_attributes(x, y, z, 2) = k_param_a * (exp (k_param_b &
                                           * layer_mass(x, y, z, 3)))
            END DO
          END IF
        END DO
      END DO Layer_k

      !Calculate transmissivity profile for each column
      CALL trans_height(x_extent, y_extent, no_layers, layer_attributes, &
                        transmissivity, activation_status)

      !4.Recalculate layer wet proportion
      Layer_wet_prop: DO x = 1, x_extent
        DO y = 1, y_extent
          IF (activation_status(x, y) == "on") THEN
            DO z = 1, (no_layers - 1)
              !Wet proportion
              IF (z == 1) THEN
                IF (water_table(x, y) >= transmissivity(x, y, z, 1)) THEN
                    wet_proportion(x, y, 1) = 1.0
                  ELSE
                    wet_proportion(x, y, 1) =  water_table(x, y) /&
                      layer_attributes(x, y, 1, 1)
                  END IF
                ELSE IF (water_table(x, y) >= transmissivity(x, y, z, 1)) THEN
                  wet_proportion(x, y, z) = 1.0
                ELSE IF (water_table(x, y) <= transmissivity(x, y, z - 1, 1)) &
                  THEN
                    wet_proportion(x, y, z) = 0.0
                ELSE
                  wet_proportion(x, y, z) = (water_table(x, y) -&
                  transmissivity(x, y, z - 1, 1)) / layer_attributes(x, y, z, 1)
              END IF
            END DO
          END IF
        END DO
      END DO Layer_wet_prop
      
    IF (week_counter == 52) EXIT

    !End sub-annual loop
    END DO Sub_annual_loop

    !8.Calculate the average water-table depth for the year
    !to be used for the next layer litter calculation
    Water_table_depth_av: DO x = 1, x_extent
      DO y = 1, y_extent
        IF (activation_status(x, y) == "on") THEN
          col_wt_depth_av(x, y) = col_wt_sum(x, y) / (REAL(sub_year_counter))
          col_wt_depth_av_su(x, y) = col_wt_sum_su(x, y) / & 
            (REAL(sub_year_counter_su))
        END IF
      END DO
    END DO Water_table_depth_av

    !Calculate the mean timestep in minutes
    mean_t_step = (1.0 / REAL(sub_year_counter)) * 525600
    WRITE(130, '(20F20.8)') mean_t_step
    !Reset mean timestep for the next main loop
    mean_t_step = 0.0

    !Reset month timeout flag
    week_counter = 0

    !6.Check for model write out
    output_counter = output_counter + 1
    IF (output_counter == output_interval) THEN
      !Write results to file
      DO x = 1, x_extent
        DO y = 1, y_extent
          IF(activation_status (x,y) == "off" &
             .OR. activation_status (x,y) == "diri" &
             .OR. activation_status (x,y) == "neu") THEN
            !Column height
            WRITE(070, '(20F20.8)') -999.0
            !Water-table height
            WRITE(080, '(20F20.8)') -999.0
            !Column mass per area
            WRITE(100, '(20F20.8)') -999.0
            !Column mean water-table depth
            WRITE(110, '(20F20.8)') -999.0
            !Summer mean column water-table depth
            WRITE(160, '(20F20.8)') -999.0

          ELSE
            WRITE(070, '(20F20.8)')&
              transmissivity(x, y, (no_layers - 1), 1)
            WRITE(080, '(20F20.8)') water_table(x, y)
            WRITE(100, '(20F20.8)') col_mass_per_area(x, y)
            WRITE(110, '(20F20.8)') col_wt_depth_av(x, y)
            WRITE(160, '(20F20.8)') col_wt_depth_av_su(x, y)
          END IF
        END DO
      END DO
    END IF

    !7.Annual counting
    !Update counter for annual update
    year_counter = year_counter + 1 ! Counts annual timesteps

    !Exit main annual loop if year counter exceeds model timespan
    IF (year_counter > t_extent) EXIT

    !Reset column water-table sum for next annual loop
    col_wt_sum = 0.0
    col_wt_sum_su = 0.0
    !Reset sub-year counting for the next annual loop
    sub_year_counter = 0
    sub_year_counter_su = 0

    !9.Add new layer for each additional year
    !current layer counter +1 for pond
    no_layers = no_layers + 1

    !Initialise above-ground storage layer properties for new ponding layer
    !for active columns only
    Pond_props: DO x = 1, x_extent
      DO y = 1, y_extent
        IF (activation_status(x, y) == "on") THEN
         !Initialise the properties of the pond
          layer_attributes(x, y, no_layers, 1) = pond_depth
          layer_attributes(x, y, no_layers, 2) = 0.0 !k
          layer_attributes(x, y, no_layers, 3) = 1.0 !s
        END IF
      END DO
    END DO Pond_props

    !10.Calculate mass of new annual litter layer using modification of
    !empirical quartic function in Belyea and Clymo (2001) Proceedings
    !Royal Society London B, 268, 1315-1321. The modification takes
    !account of temperature effects.

    !Is average water table > 66.8 cm below surface? If so, assign v. low
    !productivity to new annual layer (= mass addition because timesteps
    !for new layer addition are annual).
    New_layer_props: DO x = 1, x_extent
      DO y = 1, y_extent
        IF (activation_status(x, y) == "on") THEN
          !Layer mass of new top layer
          IF (col_wt_depth_av(x, y) > 66.8) THEN
            layer_mass(x, y, (no_layers - 1), 1) = 0.0000001
          ELSE
            layer_mass(x, y, (no_layers - 1), 1) =&
                        ((9.3 + (1.33 * col_wt_depth_av(x, y)) - (0.022&
                        * (col_wt_depth_av(x, y)**2.0)))**2.0) * 0.0001&
                        !Use mean annual value from monthly weather
                        *((0.1575 * mean_temp) + 0.009132)
          END IF
            !Layer initial mass of new top layer = layer mass of top layer
            layer_mass(x, y, (no_layers - 1), 2) =&
              layer_mass(x, y, (no_layers - 1), 1)
            !Layer remaining mass of new top layer = 100%
            layer_mass(x, y, (no_layers - 1), 3) = 1.0
            !Layer thickness of new top layer
            layer_attributes(x, y, (no_layers - 1), 1) =&
              layer_mass(x, y, (no_layers - 1), 1) / density
            !Layer drainable porosity
            layer_attributes(x, y, (no_layers - 1), 3) = porosity
            !Layer elevation of new top layer =
            transmissivity(x, y, (no_layers - 1), 1) =&
              !layer elevation below new top layer + thickness of new layer
              transmissivity(x, y, (no_layers - 1) - 1, 1) +&
                layer_attributes(x, y, (no_layers - 1), 1)
            !Wet proportion of new top layer
            IF (water_table(x, y) >= transmissivity(x, y, (no_layers - 1), 1)) &
              THEN
                wet_proportion(x, y, (no_layers - 1)) = 1.0
              ELSE IF (water_table(x, y) <= &
                transmissivity(x, y, (no_layers - 1), 1)) THEN
                  wet_proportion(x, y, z) = 0.0
              ELSE
                wet_proportion(x, y, (no_layers -1)) = (water_table(x, y) -&
                  transmissivity(x, y, (no_layers - 1), 1)) &
                  / layer_attributes(x, y, (no_layers - 1), 1)
            END IF
            !Set layer Ksat if new layer contains water 
            IF (wet_proportion(x, y, (no_layers - 1)) > 0.0) THEN
              layer_attributes(x, y, (no_layers - 1), 2) =  &
                k_param_a * (exp (k_param_b))
            !Calculate transmissivity with new layer
              transmissivity(x, y, (no_layers - 1), 2) = &
                transmissivity(x, y, (no_layers - 1) - 1, 2) & 
                + (layer_attributes(x, y, (no_layers -1), 1) &
                * layer_attributes(x, y, (no_layers - 1), 2))
            END IF
            col_mass_per_area(x, y) = col_mass_per_area(x, y) +&
              layer_mass(x, y, (no_layers - 1), 1)
        END IF
      END DO
    END DO New_layer_props

   !Re-set output counter
   output_counter = 0

  !End of main annual loop
  END DO Annual_loop

  !11.Write array outputs
  DO x = 1, x_extent
    DO y = 1, y_extent
      IF (activation_status(x, y) == "on") THEN
        DO z = 2, (no_layers -1)
          WRITE(060, '(20F20.8)') layer_mass(x, y, z, 3)     !mass remaining
          WRITE(090, '(20F20.8)') transmissivity(x, y, z, 2) !col transmiss'y
          WRITE(120, '(20F20.8)') wet_proportion(x, y, z)    !layer wet prop_n
          WRITE(150, '(20F20.8)') layer_attributes(x, y, z, 2) !layer k
        END DO
      END IF
    END DO
  END DO

END PROGRAM DIGIBOG_BB
