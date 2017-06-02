# TsunamiGlobe
C++ library for modeling tsunamis from Virtual Quake simulated earthquakes using the Tsunami Squares method (Steven N. Ward).  Extended from c++ code written by Kasey W. Schultz and John Max Wilson with intent to model a spherical or ellipsoidal Earth and improve stability.  NOAA ETOPO1 combined bathymetry and topology data included.
 

To compile, use
$ bash setup.sh

Run simulation with
$ ./TsunamiSquares tsunami_params.txt

Visualizations can be created with
$ python tsunamisquares.py


== Parameter Definitions ==
tsunami_params.txt contains simulation parameters which may be changed without recompiling from source.

Parameter						Definition
--------------------------------------------------
out\_file\_name 					File name for simulation result output

bathy_file 						File with bathymetry/topography data

kml_file 						File for KML output (currently unused)

deformation_file 				File name for initial seafloor uplift data (e.g. from thrust earthquake)

move_bool						Whether to do the main square-moving simulation

accel_bool						Whether to accellerate water each step

diffuse_bool					Whether to diffuse (smooth) water each step

D 								Diffusion constant

dt								Time step in seconds

N_steps 						Number of time steps in simulation

current_step 					Initial time step

update_step 					Period of printing simulation state to terminal

save_step 						Period of printing simulation state to file

time 							Initial time

output\_num\_digits\_for\_percent 	Float precision for output

flatten_bool					Whether to flatten the seafloor before running, for testing purposes

flat_depth 						Used when flattening the bathymetry to a constant depth for testing purposes (negative for below sea level)

bump_bool						Whether to use bump instead deformation file, for testing purposes

bump_height 					Used when testing, produces small initial seafloor displacement
