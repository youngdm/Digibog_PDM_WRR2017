MODULE hydro_procedures_bb_transmax

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


!-------------------------------------------------------------------------------
! Section 1.0 Program header
!-------------------------------------------------------------------------------

  !For details of code ownership, model updates etc, see section 1.0 of
  !DigiBog_Hydro
  
  !-----------------------------------------------------------------------------
  !Modification history of code 
  !-----------------------------------------------------------------------------
  !Programmer           Date           Modifications
  !-----------------------------------------------------------------------------
  !Dylan M. Young       18/06/14  This version is for the fixed Dirichlet      
  !                               condition model. The activation status in  
  !                               trans_height (section 2.0) and wat_k_mean
  !                               (section 4.0) have been removed because the 
  !                               'diri' condition is updated using only water-
  !                               table height.
  !-----------------------------------------------------------------------------
  !Dylan M. Young       25/06/14  Transmax scalar added to wat_K_mean to enable
  !                               variable timesteps to be used. Transmax code 
  !                               created by Andy Baird.  
  !-----------------------------------------------------------------------------
  !Dylan M. Young       23/01/15  Water-table-update (6.0). Added condition to
  !                               set marker = no_layers where the water-table
  !                               exceeds the ponding layer. The layer marker
  !                               is initialised at 0.
  !-----------------------------------------------------------------------------
  !Declarations

  IMPLICIT NONE

  !Global type definition
  !Define a real kind type q with at least 8 decimal digits and an exponent
  !range from 10**30 to 10**(-30)
  INTEGER, PARAMETER :: q = SELECTED_REAL_KIND(P = 8, R = 30)

  CONTAINS


!-------------------------------------------------------------------------------
! Section 2.0   Initial transmissivity calculations
!-------------------------------------------------------------------------------

  SUBROUTINE trans_height(x_extent, y_extent, no_layers, layer_attributes, &
                          transmissivity, activation_status)

    !This subroutine calculates the peat column's transmissivity from the
    !impermeable base of each column to the top of each layer


    !Declarations

    IMPLICIT NONE 

    !Subroutine arguments 

    !Scalar arguments with INTENT(IN) 
    INTEGER, INTENT(IN) :: x_extent, y_extent, no_layers
 
    !Array arguments with INTENT(IN)
    REAL(KIND=q), DIMENSION(:,:,:,:), INTENT(IN) :: layer_attributes
    CHARACTER(8), DIMENSION(:,:), INTENT(IN) :: activation_status

    !Array arguments with INTENT(INOUT)
    REAL(KIND=q), DIMENSION(:,:,:,:), INTENT(INOUT) :: transmissivity

    !Local scalars
    INTEGER :: x, y, z 
         
    !-- End of subroutine header -----------------------------------------------

    !Calculations
    
    DO x = 1, x_extent
      DO y = 1, y_extent
        !Calculations only needed for active columns and Dirichlet boundary
        !columns
        IF(activation_status(x, y)     == "on") THEN 
          transmissivity(x, y ,1, 1) = layer_attributes(x, y, 1, 1)
          transmissivity(x, y, 1, 2) = layer_attributes(x, y, 1, 1) &
                                     * layer_attributes(x, y, 1, 2)        
          DO z = 2, (no_layers - 1)!Ignore ponding layer
            transmissivity(x, y, z, 1) = transmissivity(x, y, z - 1, 1) &
                                       + layer_attributes(x, y, z, 1)
            transmissivity(x, y, z, 2) = transmissivity(x, y, z - 1, 2) & 
                                       + (layer_attributes(x, y, z, 1) &
                                       * layer_attributes(x, y, z, 2))
          END DO
        END IF
      END DO
    END DO

  END SUBROUTINE trans_height


!-------------------------------------------------------------------------------
! Section 3.0   Water capacity initialisation
!-------------------------------------------------------------------------------

  SUBROUTINE layer_water_depth (x_extent, y_extent, no_layers, &
                                layer_attributes, layer_storage, &
                                activation_status)

    !This subroutine calculates the amount of water (expressed as a depth) 
    !capable of being stored in each layer in each column.

    !Declarations

    IMPLICIT NONE 

    !Subroutine arguments
    
    !Scalar arguments with INTENT(IN) 
    INTEGER, INTENT(IN) :: x_extent, y_extent, no_layers

    !Array arguments with INTENT(IN) 
    REAL(KIND=q), DIMENSION(:,:,:,:), INTENT(IN) :: layer_attributes
    CHARACTER(8), DIMENSION(:,:), INTENT(IN) :: activation_status
  
    !Array arguments with INTENT(INOUT) 
    REAL(KIND=q), DIMENSION(:,:,:,:), INTENT(INOUT) :: layer_storage

    !Local scalars
    INTEGER :: x, y, z
   
    !-- End of subroutine header -----------------------------------------------

    !Calculations
    
    DO x = 1, x_extent
      DO y = 1, y_extent
        !Calculations only needed for active columns
        IF (activation_status(x, y) == "on") THEN
          layer_storage(x, y, 1, 1) = layer_attributes(x, y, 1, 1)
          layer_storage(x, y, 1, 2) = layer_attributes(x, y, 1, 1) &
                                    * layer_attributes(x, y, 1, 3)        
          DO z = 2, no_layers !Ponding layer included
            layer_storage(x, y, z, 1) = layer_storage(x, y, z - 1, 1) & 
                                      + layer_attributes(x, y, z ,1)
            layer_storage(x, y, z, 2) = layer_attributes(x, y, z, 1) &
                                      * layer_attributes(x, y, z, 3)
          END DO
        END IF
      END DO
    END DO

  END SUBROUTINE layer_water_depth


!------------------------------------------------------------------------------
!Section 4.0 Depth-averaged hydraulic conductivity calculations
!------------------------------------------------------------------------------

  SUBROUTINE wat_k_mean(x_extent, y_extent, no_layers, transmax, water_table, &
                      wk_mean, layer_attributes, transmissivity, &
                      activation_status)

  !This subroutine calculates the depth-averaged K below the water-table
  !for each column. It also records the value of the transmissivity for
  !that model column with the highest transmissivity below the
  !water-table for use in the variable time-step calculation in the main
  !programme.

  !Declarations                    

  IMPLICIT NONE

  !Subroutine arguments
  
  !Scalar arguments with INTENT(IN)
  INTEGER, INTENT(IN) :: x_extent, y_extent, no_layers
                    
  !Array arguments with INTENT(IN)
  REAL(KIND=q), DIMENSION(:,:), INTENT(IN) :: water_table
  REAL(KIND=q), DIMENSION(:,:,:,:), INTENT(IN) :: layer_attributes, &
                                                  transmissivity
  CHARACTER(8), DIMENSION(:,:), INTENT(IN) :: activation_status
  
  !Scalar arguments with INTENT(INOUT)
  REAL(KIND=q), INTENT(INOUT) :: transmax

  !Array arguments with INTENT(INOUT)
  REAL(KIND=q), DIMENSION(:,:), INTENT(INOUT) :: wk_mean

  !Local scalars
  INTEGER :: x, y, z

  !-- End of subroutine header-----------------------------------------------
  
  !Calculations
  DO x = 1, x_extent
    DO y= 1, y_extent
      !Calculations only needed for active columns
      IF (activation_status(x, y) == "on") THEN
        !If the water-table is at or above the peatland surface then mean K is
        !that of the peat profile below the ponding layer
        IF (water_table(x, y) &
            >= transmissivity(x, y, (no_layers -1), 1)) THEN
          wk_mean(x, y) = transmissivity(x, y, (no_layers -1), 2)
          IF (wk_mean(x, y) > transmax) THEN
            transmax = wk_mean(x, y)
          END IF
          wk_mean(x, y) = wk_mean(x, y) / & 
                            transmissivity(x, y, (no_layers -1), 1)
        !If the water-table is below the top of the first layer
        ELSE IF (water_table(x, y) < transmissivity(x, y, 1 , 1)) THEN
          wk_mean(x, y) = layer_attributes(x, y, 1, 2)
          IF (transmissivity(x, y, 1, 2) > transmax) THEN
            transmax = transmissivity(x, y, 1, 2)
          END IF
        !If the water-table is equal to or above the top of the first layer and
        !below the top of the bog surface
        ELSE
          DO z = 1, (no_layers - 2)
            !If the water-table is equal to the height of layer z
            IF (water_table(x, y) == transmissivity(x, y, z, 1)) THEN
              wk_mean(x, y) = transmissivity(x, y, z, 2)
              IF (wk_mean(x, y) > transmax) THEN
                transmax = wk_mean(x, y)
              END IF
              wk_mean(x, y) = wk_mean(x, y) / transmissivity(x, y, z, 1)
            !If the water-table is above the top of layer z
            ELSE IF (water_table(x, y) > transmissivity(x, y, z, 1)) THEN
            !If the water-table is below the top of z + 1
              IF (water_table(x, y) < transmissivity(x, y, z + 1, 1)) THEN
              wk_mean(x, y) = (transmissivity(x, y, z, 2) &
                            + (layer_attributes(x, y, z + 1, 2) &
                            * (water_table(x, y) &
                            - transmissivity(x, y, z, 1))))
                IF (wk_mean(x, y) > transmax) THEN
                transmax = wk_mean(x, y)
                END IF
              wk_mean(x, y) = wk_mean(x, y) / water_table(x, y)
              END IF
            END IF
          END DO
        END IF
      END IF
    END DO
  END DO

  END SUBROUTINE wat_k_mean


!-------------------------------------------------------------------------------
! Section 5.0   2-dimensional flux calculations
!-------------------------------------------------------------------------------

  SUBROUTINE move_water(x_extent, y_extent, timestep, spatial_step, rainfall, &
                        base_altitude, water_change, water_table, wk_mean, &
                        activation_status)
  
    !This subroutine calculates the net movement of water (expressed as a depth) 
    !between columns. It is assumed that a harmonic mean operates between the
    !mean (depth-averaged) K of each column. The output is a change in depth
    !(positive or negative) of water (volume per unit area)to be added to the
    !previous water-table elevation.
    
    !Declarations

    IMPLICIT NONE   
                        
    !Subroutine arguments
    
    !Scalar arguments with INTENT(IN)
    INTEGER, INTENT(IN) :: x_extent, y_extent
    REAL(KIND=q), INTENT(IN) :: spatial_step, timestep, rainfall
  
    !Array arguments with INTENT(IN) 
    REAL(KIND=q), DIMENSION(:,:), INTENT(IN) :: base_altitude, water_table, &
                                                wk_mean
    CHARACTER(8), DIMENSION(:,:), INTENT(IN) :: activation_status

    !Array arguments with INTENT(INOUT)                                               
    REAL(KIND=q), DIMENSION(:,:), INTENT(INOUT) :: water_change

    !Local scalars
    INTEGER :: x, y   

    !Local arrays 
    REAL(KIND=q), DIMENSION(x_extent, y_extent) :: x_flux, y_flux

                         
    !-- End of subroutine header -----------------------------------------------


    !Calculations
    
    !Volume 'flux' (in x direction) between columns using Dupuit-Forchheimer
    !approximation
    !16.09.2013 Changed way in which flow to and from Diri cells works
    DO x = 1, (x_extent - 1)
      DO y = 2, (y_extent - 1)
        !Different rules apply to cells with different statuses
        !Case 1
        IF (activation_status(x, y) == "neu" &
           .AND. activation_status(x + 1, y) == "on") THEN 
          x_flux(x, y) = 0.0
        !Case 2
        ELSE IF (activation_status(x, y) == "diri" &
                 .AND. activation_status(x + 1, y) == "on") THEN
          x_flux(x, y) = wk_mean(x + 1, y) &
                       * ((water_table(x, y) + water_table(x + 1, y)) / 2.0) &
                       * (base_altitude(x, y) + water_table(x, y) &
                       - base_altitude(x + 1, y)  - water_table(x + 1, y)) &
                       * 2.0 * timestep
        !Case 3
        ELSE IF (activation_status(x, y) == "on" &
                 .AND. activation_status(x + 1, y) == "neu") THEN
          x_flux(x, y) = 0.0
        !Case 4
        ELSE IF (activation_status(x, y) == "on" &
                 .AND. activation_status(x + 1, y) == "diri") THEN
          x_flux(x, y) = wk_mean(x, y) &
                       * ((water_table(x, y) + water_table(x + 1, y)) / 2.0) &
                       * (base_altitude(x, y) + water_table(x, y) &
                       - base_altitude(x + 1, y)  - water_table(x + 1, y)) &
                       * 2.0 * timestep
        !Case 5
        ELSE IF (activation_status(x, y) == "on" &
                 .AND. activation_status(x + 1, y) == "on") THEN
          x_flux(x, y) = (2 * wk_mean(x, y) * wk_mean(x + 1, y)) &
                       / (wk_mean(x, y) + wk_mean(x + 1, y)) &
                       * ((water_table(x, y) + water_table(x + 1, y)) / 2.0) &
                       * (base_altitude(x, y) + water_table(x, y) &
                       - base_altitude(x + 1, y)  - water_table(x + 1, y)) &
                       * timestep
        END IF
      END DO
    END DO

    !Volume 'flux' (in y direction) between columns using Dupuit-Forchheimer
    !approximation
    DO y = 1, (y_extent - 1)
      DO x = 2, (x_extent - 1)
        !Different rules apply to cells with different statuses
        !Case 1
        IF (activation_status(x, y) == "neu" &
           .AND. activation_status(x, y + 1) == "on") THEN 
          y_flux(x, y) = 0.0
        !Case 2
        ELSE IF (activation_status(x, y) == "diri" &
                 .AND. activation_status(x, y + 1) == "on") THEN
          y_flux(x, y) = wk_mean(x, y + 1) &
                       * ((water_table(x, y) + water_table(x, y + 1)) / 2.0) &
                       * (base_altitude(x, y) + water_table(x, y) &
                       - base_altitude(x, y + 1)  - water_table(x, y + 1)) &
                       * 2.0 * timestep
        !Case 3
        ELSE IF (activation_status(x, y) == "on" &
                 .AND. activation_status(x, y + 1) == "neu") THEN
          y_flux(x, y) = 0.0
        !Case 4
        ELSE IF (activation_status(x, y) == "on" &
                 .AND. activation_status(x, y + 1) == "diri") THEN
          y_flux(x, y) = wk_mean(x, y) &
                       * ((water_table(x, y) + water_table(x, y + 1)) / 2.0) &
                       * (base_altitude(x, y) + water_table(x, y) &
                       - base_altitude(x, y + 1)  - water_table(x, y + 1)) &
                       * 2.0 * timestep
        !Case 5
        ELSE IF (activation_status(x, y) == "on" &
                 .AND. activation_status(x, y + 1) == "on") THEN
          y_flux(x, y) = (2 * wk_mean(x, y) * wk_mean(x, y + 1)) &
                       / (wk_mean(x, y) + wk_mean(x, y + 1)) &
                       * ((water_table(x, y) + water_table(x, y + 1)) / 2.0) &
                       * (base_altitude(x, y) + water_table(x, y) &
                       - base_altitude(x, y + 1)  - water_table(x, y + 1)) &
                       * timestep
        END IF
      END DO
    END DO

    !Convert volume into depth of water (in absence of soil matrix - i.e. the
    !depth of water if the drainable porosity were 1)
    !Ignore all edge cells (they have to be boundary conditions)
    DO x = 2, (x_extent - 1)
      DO y = 2, (y_extent - 1)
        IF (activation_status(x, y) == "on") THEN
          water_change(x, y) = (x_flux(x - 1, y) - x_flux(x ,y) &
                             +  y_flux(x, y - 1) - y_flux(x ,y)) &
                             / (spatial_step ** 2)
          water_change(x, y) = water_change(x, y) + (rainfall * timestep)
        END IF
      END DO
    END DO

  END SUBROUTINE move_water


!-------------------------------------------------------------------------------
! Section 6.0   Water-table update
!-------------------------------------------------------------------------------

  SUBROUTINE water_table_update(x_extent, y_extent, no_layers, water_change, &
                                water_table, layer_attributes, layer_storage, &
                                activation_status) 
 
    !This subroutine updates the position of the water table based on the 
    !storage available in layers above or below the water table
    
    !Declarations

    IMPLICIT NONE  

    !Subroutine arguments
    
    !Scalar arguments with INTENT(IN)
    INTEGER, INTENT(IN) :: x_extent, y_extent, no_layers 
         
    !Array arguments with INTENT(IN) 
    REAL(KIND=q), DIMENSION(:,:,:,:), INTENT(IN) :: layer_attributes, &
                                                    layer_storage
    CHARACTER(8), DIMENSION(:,:), INTENT(IN) :: activation_status
 
    !Array arguments with INTENT(INOUT):
    REAL(KIND=q), DIMENSION(:,:), INTENT(INOUT) :: water_change, water_table
         
    !Local scalars 
    INTEGER :: x, y, z, marker
    REAL(KIND=q) :: test_depth
  
    !-- End of subroutine header -----------------------------------------------
                                 
    !Calculations
    
    !Ignore all edge cells (they have to be boundary conditions)
    DO x = 2, (x_extent - 1)
      DO y = 2, (y_extent - 1)
        IF (activation_status(x, y) == "on") THEN
          !Locate layer in which water table resides
          IF (water_table(x, y) <= layer_storage(x, y, 1, 1)) THEN
            marker = 1
          ELSE
            DO z = 1, (no_layers - 1)
              IF (water_table(x, y) > layer_storage(x, y, z, 1)) THEN
                !If water table is in layer z + 1
                IF (water_table(x, y) <= layer_storage(x, y, z + 1, 1)) THEN
                  marker = z + 1
                ELSE IF (water_table(x, y) > &
                  layer_storage(x, y, no_layers, 1)) THEN
                    marker = no_layers
                  EXIT
                END IF
              END IF
            END DO
          END IF
          !Does water table rise?
          IF (water_change(x, y) > 0.0) THEN 
            !Ignore any water-table updates if water level already at top
            !of uppermost cell
            IF (water_table(x, y) < layer_storage(x, y, no_layers, 1)) THEN 
              z = marker
              !If water table at top of layer z 
              IF (water_table(x, y) == layer_storage(x, y, z, 1)) THEN
                !No storage available - advance to next layer
                z = z + 1
              END IF
              DO
                test_depth = layer_storage(x, y, z, 2) &
                           * (layer_storage(x, y, z, 1) - water_table(x, y)) &
                           / layer_attributes(x, y, z, 1)
                IF (water_change(x, y) <= test_depth) THEN
                  !Water table moves up within layer z
                  water_table(x, y) = water_table(x, y) + (water_change(x, y) &
                                    / layer_attributes(x, y, z, 3))
                  EXIT
                !Water table moves up to next layer or rises to the surface
                ELSE
                  !Has top of ponding layer been reached?
                  IF (z == no_layers) THEN
                    water_table(x, y) = layer_storage(x, y, z, 1)
                    EXIT
                  ELSE
                    water_table(x, y) = layer_storage(x, y, z, 1)
                    water_change(x, y) = water_change(x, y) - test_depth
                    z = z + 1
                  END IF
                END IF              
              END DO
            END IF
          END IF
          !Does the water table fall?
          IF (water_change(x, y) < 0.0) THEN
            !Ignore any water-table updates if water table already at or 
            !below base of column
            IF (water_table(x, y) > 0.0 ) THEN
              z = marker
              DO 
                IF (z == 1) THEN
                  test_depth = layer_storage(x, y, 1, 2) * water_table(x ,y) &
                             / layer_attributes(x, y, 1, 1)
                ELSE
                  test_depth = layer_storage(x, y, z, 2) & 
                             * (water_table(x, y) &
                             - layer_storage(x, y, z - 1, 1)) &
                             / layer_attributes(x, y, z, 1)
                END IF
                IF (ABS(water_change(x, y)) <= test_depth) THEN 
                  !Water table falls within layer z or to top of layer z - 1
                  water_table(x, y) = water_table(x, y) + (water_change(x, y) &
                                    / layer_attributes(x, y, z, 3) )
                  EXIT
                ELSE !Water table falls to top of layer z - 1 (if it exists)
                     !and needs to fall further, provided a lower layer exits
                  !Has base of column been reached/exceeded?
                  IF (z == 1) THEN
                    water_table(x, y) = 0.0
                    !If this were to happen might show model instabilty, so
                    !perhaps add code to stop model run
                    EXIT
                  ELSE
                    water_table(x, y) = layer_storage(x, y, z - 1, 1)
                    water_change(x, y) = water_change(x, y) + test_depth
                    z = z - 1
                  END IF
                END IF              
              END DO
            END IF
          END IF
        END IF
      END DO
    END DO

  END SUBROUTINE water_table_update
  
END MODULE hydro_procedures_bb_transmax
