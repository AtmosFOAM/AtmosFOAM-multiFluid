#!/bin/bash

# Create data files for single column, multi-fluid cases, consistent with
# horizontal means of other cases

if [ "$#" -ne 2 ]
then
   echo usage: ./plotProfiles.sh caseName time
   exit
fi

case=$1/hMean
time=$2

# A bit more file manipulation, Write out ascii data and sort by z
for part in stable buoyant; do
    if [ -e $case/$time/u.$part ]; then
        # Write out components of velocity and rename
        writeuvw u.$part -case $case -time $time
        mv $case/$time/u.${part}z $case/$time/uz.${part}
        # Calculate P.part
        sumFields -case $case $time P.$part $time P $time Pi.$part
    fi

    for var in b uz P sigma massTransfer.buoyant massTransfer.stable; do
        if [ -e $case/$time/$var.$part ]; then
            writeCellDataxyz -case $case -time $time $var.$part
                sort -g -k 3 $case/$time/$var.$part.xyz \
                    | sponge $case/$time/$var.$part.xyz
        fi
    done
done
writeCellDataxyz -case $case -time $time b
    sort -g -k 3 $case/$time/b.xyz \
        | sponge $case/$time/b.xyz

# Plots
for var in b P S sigma w; do
    (cd $case/$time && pwd && gmtPlot ../../plots/$var.gmt)
done


