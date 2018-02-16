#!/bin/bash

g++ -I /usr/local/include/boost_1_64_0 main.cpp TsunamiObjects.cpp TsunamiUtil.cpp -o TsunamiGlobe -lGeographic -lnetcdf_c++4 -lalglib -fopenmp
