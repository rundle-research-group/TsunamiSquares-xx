#!/bin/bash

g++ -I /usr/local/include/boost_1_64_0 main.cpp TsunamiGlobe.cpp TsunamiGlobeUtil.cpp -o TsunamiGlobe -l Geographic -fopenmp
