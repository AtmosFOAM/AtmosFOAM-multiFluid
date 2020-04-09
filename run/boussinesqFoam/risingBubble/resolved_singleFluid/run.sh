#!/bin/bash -e

case=.

# clear out old stuff
rm -rf $case/[0-9]* $case/constant/polyMesh $case/core $case/log

# create mesh
blockMesh -case $case

# Stationary initial conditions
mkdir $case/0
cp -r $case/init_0/* $case/0

# Add a buoyant perturbation
setInitialTracerField -name b -tracerDict setbDict -case $case

# Solve Boussinesq equations (boussinesqFoam also writes out the initial
# velocity plua a tiny increment for conditional averaging)
boussinesqFoam -case $case >& log & sleep 0.01; tail -f log

# Plot the final result
time=1000
writeuvw -time $time u -case $case
gmtFoam -time $time bU -case $case
gv $case/$time/bU.pdf &
gmtFoam -time $time bP -case $case
gv $case/$time/bP.pdf &

# Conditional average various time steps
rm -r $case/hMean/[0-9]*
blockMesh -case $case/hMean
for time in [0-9]*; do
    ./conditionalAverage.sh $case $time
done

# Plots
for time in 0 500 1000; do
    plots/plotProfiles.sh $case $time
done

