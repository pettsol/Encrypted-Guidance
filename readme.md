# Encrypted Guidance

## Overview

This repository contains source code and instructions for a proof-of-concept simulation for an encrypted guidance system. The encrypted guidance system computes encrypted yaw references using a linearized LOS guidance law without knowing of the vehicles' position and desired destination.

** Author: [Petter Solnoer](https://www.ntnu.no/ansatte/petter.solnor), petter.solnor@ntnu.no **

The simulation can be compiled with g++ and executed in Linux with the following commands:

Compile: g++ main.cpp he\_guidance.cpp encoder.cpp joye\_libert\_journal/joye\_libert.cpp -o main -lgmp -lgmpxx
Execute: ./main

The implementation uses the GMP big number library. For more information, please visit https://gmplib.org/
