#!/bin/bash
profilesDistance="1 2"
spotDistance="1 2"
spotDistanceTransformation="none linear log sqrt"
for p in $profilesDistance
do 
    for sd in $spotDistance
    do
        for sdt in $spotDistanceTransformation
        do
            Rscript --vanilla ./cluster.R workingDir="../Datasets/HumanLymphNode" \
                profilesDistance=$p spotDistance=$sd spotDistanceTransformation=$sdt
        done
    done
done