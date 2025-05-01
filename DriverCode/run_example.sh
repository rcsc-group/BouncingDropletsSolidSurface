#!/bin/bash

# Author: Radu Cimpeanu
# Date: 24/10/2024

# Additional resolution levels, drop radii or velocities can be added below

for LEVEL in 12; do
	for RADIUS in 0.00021; do
		for VELOCITY in 0.481913; do

			cp -r MasterImpact/ ContactRadius-R$RADIUS-V$VELOCITY-Level$LEVEL
			cd ContactRadius-R$RADIUS-V$VELOCITY-Level$LEVEL/

			qcc -O2 -w -fopenmp -Wall DropImpact.c -lm -o DropImpact -L$BASILISK/gl -lglutils -lfb_tiny
			# parameters:
			# 1. rho liquid
			# 2. rho gas
			# 3. mu liquid
			# 4. mu gas
			# 5. sigma (surface tension coefficient)
			# 6. g acceleration
			# 7. drop radius
			# 8. initial drop velocity
			# 9. simulation end time
			# 10. max level

			export OMP_NUM_THREADS=4

			# send the three parameters in the sweep, the level and max runtime
			# For this case, no gravity and 
			./DropImpact 870 1.21 0.00174 1.81e-5 0.0187 9.81 $RADIUS $VELOCITY 10.001 $LEVEL

			cd ..
		
		done
	done
done
