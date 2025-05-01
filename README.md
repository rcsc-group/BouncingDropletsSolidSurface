# Bouncing droplets onto a solid substrate
Direct numerical simulation code infrastructure for drop impact onto solid substrates in the inertio-capillary regime, supporting collaborative work with the Harris Lab at Brown. This represents the sister implementation to the impact onto liquid substrate scenario which can be found under the [BouncingDroplets](https://github.com/rcsc-group/BouncingDroplets/) repository. 

## Installation
* The code relies on [Basilisk](<http://basilisk.fr/>) to model the Navier-Stokes equations. See the [installation page](<http://basilisk.fr/src/INSTALL>) for instructions. 
* Full visualisation capabilities have been used in order to generate animations. These may be switched off depending on the local architecture.

## Running the code
Once the Basilisk structure is in place, the driver code here is built in order to navigate parameter sweeps in velocity $V_0$ and resolution level, with one of each values added to the run_master_example.sh for brevity. Other parameters can be varied through this shell script, with both physical and computational handles provided. 

The code can be executed by simply executing this shell script via *sh run_master_example.sh* inside a terminal. Output will then be produced within a foldering structure that consists of summary DNS execution information, mass conservation and VOF data, interface coordinates, simulation slices and animations, which can be used for further post-processing.
