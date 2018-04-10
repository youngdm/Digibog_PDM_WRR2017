# Digibog_PDM_WRR2017
DigiBog peatland development model code used in the paper 'Simulating the long‐term impacts of drainage and restoration on the ecohydrology of peatlands' https://doi.org/10.1002/2016WR019898. The model was used to simulate the growth and subsequent drainage of a sloping peat transect. 

This version of the model runs in 2.5-D and comprises two files written in Fortran 90. They are;
DigiBog_WRR2017_main.f90 and DigiBog_WRR2017_hydro_module.f90. A drain is not simulated in the main model code included here, but can be added to an active column by changing its status to a Dirichlet boundary condition at some point in the model run, and then fixing a water-table position to the required drain depth. The damming of the drain can be simulated by raising the water-table in the drained column. 

The model uses five input files which are also included with the model code;
1. User-defined input parameters - 010_DigiBog_BB_IN_information.txt
2. Net rainfall time series* - 020_DigiBog_BB_IN_net_rain.txt
3. Temperature time series* - 030_DigiBog_BB_IN_temp.txt
4. Column status (active, boundary column, or inactive) - 040_DigiBog_BB_IN_column_status.txt
5. Base altitude for peat columns - 050_DigiBog_BB_IN_baltitude.txt

*Note these net rainfall and temperature file time series have been randomly generated and are not the same files that were used in the publication.
