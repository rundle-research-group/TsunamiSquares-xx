# TsunamiGlobe
C++ library for modeling tsunamis from Virtual Quake simulated earthquakes using the Tsunami Squares method (Steven N. Ward).  Extended from c++ code written by Kasey W. Schultz and John Max Wilson (<https://github.com/johnmaxwilson/TsunamiSquares>) with intent to model a spherical or ellipsoidal Earth and improve stability.

## Requirements
NOAA ETOPO1 combined bathymetry and topology data not included, can be downloaded from <https://www.ngdc.noaa.gov/mgg/global/>  
Generating seafloor uplift initial conditions requires the quakelib library, part of Virtual Quake (<https://geodynamics.org/cig/software/vq/>).

Additionally, the following packages and libraries are required
### C++
OpenMP  
Geographiclib  
boost v1.64  
netCDF  
ALGLIB  

### Python
python 2.7  
numpy  
matplotlib  
Basemap  
ffmpeg  
scipy  
Geographiclib  
netCDF4  

### How to Use
To compile, do  
`$ bash setup.sh`

Run simulation with  
`$ ./TsunamiGlobe tsunami_params.txt`

Visualizations, bathymetry subsets, and initial condition displacement fields can be created with tsunami_tools.py.  For usage options, do  
`$ python tsunami_tools.py --help`


## Parameter Definitions  
tsunami_params.txt contains simulation parameters which may be changed without recompiling from source.

|  Parameter | Definition |
|:-----------|-----------|
|out\_file\_name | File name for simulation result output|
|bathy\_file | File with bathymetry / topography data|
|kml\_file | File for KML output (currently unused)|
|initial\_conditions | Type of initial conditions for simulation. Choose from eq, bump, gauss, and saved|
|deformation\_file | File name for initial seafloor uplift data from simulated earthquake (use with initial\_conditions  = eq)|
|bump\_height | Used when testing, produces small initial seafloor displacement (use with initial\_conditions  = bump)|
|gauss\_height | Used when testing, height of initial gaussian pile of water (use with initial\_conditions  = gauss)|
|gauss\_std | Used when testing, width of initial gaussian pile of water (use with initial\_conditions  = gauss)|
|initialstate\_file\_name | File name for initial conditions, loaded from previous simulation (use with initial\_conditions  = saved)|
|write\_sim\_state | Whether to write final simulation state at end of simulation|
|finalstate\_file\_name	| File name for final simulation state save|
|move\_bool | Whether to do the main square-moving simulation|
|accel\_bool | Whether to accellerate water each step|
|diffuse\_bool | Whether to diffuse (smooth) water each step|
|ndiffusions | How many smoothing sweeps to perform each step|
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
|doPlaneFit | Whether to use the plane-fitting method for acceleration calculations|
|num\_threads | How many threads to use during multiprocessing blocks|
|check\_sim\_health | Whether to check for NaNs and infs each time step and break the sim with helpful info if found |
