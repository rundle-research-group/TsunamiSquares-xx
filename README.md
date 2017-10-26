# TsunamiGlobe
C++ library for modeling tsunamis from Virtual Quake simulated earthquakes using the Tsunami Squares method (Steven N. Ward).  Extended from c++ code written by Kasey W. Schultz and John Max Wilson (<https://github.com/johnmaxwilson/TsunamiSquares>) with intent to model a spherical or ellipsoidal Earth and improve stability.  NOAA ETOPO1 combined bathymetry and topology data not included, can be downloaded from <https://www.ngdc.noaa.gov/mgg/global/>
 
All tools for generating initial seafloor uplift initial conditions use the quakelib library, part of Virtual Quake (<https://geodynamics.org/cig/software/vq/>).


To compile, use  
'$ bash setup.sh'

Run simulation with  
'$ ./TsunamiSquares tsunami_params.txt'

Visualizations can be created with  
'$ python tsunamisquares.py'


###Parameter Definitions  
tsunami_params.txt contains simulation parameters which may be changed without recompiling from source.

|  Parameter | Definition |  
|:-----------|-----------|  
|out\_file\_name | File name for simulation result output|
|bathy\_file | File with bathymetry / topography data|
|kml\_file | File for KML output (currently unused)|
|deformation_file | File name for initial seafloor uplift data (e.g. from thrust earthquake)|
|move\_bool | Whether to do the main square-moving simulation|
|accel\_bool | Whether to accellerate water each step|
|ndiffusions | How many smoothing sweeps to perform each step|
|diffuse\_bool | Whether to diffuse (smooth) water each step|
|D | Diffusion constant|
|dt | Time step in seconds|
|N\_steps | Number of time steps in simulation|
|current\_step | Initial time step|
|update\_step | Period of printing simulation state to terminal|
|save\_step | Period of printing simulation state to file|
|time | Initial time|
|output\_num\_digits\_for\_percent |Float precision for output|
|flatten\_bool | Whether to flatten the seafloor before running, for testing purposes|
|flat\_depth | Used when flattening the bathymetry to a constant depth for testing purposes (negative for below sea level)|
|bump\_bool | Whether to use bump instead deformation file, for testing purposes|
|bump\_height | Used when testing, produces small initial seafloor displacement|
|num\_nearest | How many neighbors to include in some nearest-neighbor overlap calculations |
|doPlaneFit | Whether to use the plane-fitting method for acceleration calculations|
|num\_threads | How many threads to use during multiprocessing blocks|
