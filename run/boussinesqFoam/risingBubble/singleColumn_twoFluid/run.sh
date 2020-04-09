#!/bin/bash -e

case=.

# clear out old stuff
rm -rf $case/[0-9]* $case/constant/polyMesh $case/core $case/log

# create mesh
blockMesh -case $case

# Initial conditions from ../resolved_singleFluid with modifications
cp -r ../resolved_singleFluid/hMean/0 .
cp 0/u 0/u.stable
cp 0/u 0/u.buoyant
sed -i 's/buoyancyf/bf.sum/g' 0/P

# Solve multi-fluid Boussinesq equations
multiFluidBoussinesqFoam >& log & sleep 0.01; tail -f log

# Plots
mv 0/P.buoyant 0/Pi.buoyant
mv 0/P.stable 0/Pi.stable
for time in 0 500 1000; do
    plots/plotProfiles.sh $case $time
done

