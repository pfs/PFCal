#!/bin/bash

siThicknesses=(120 200 320)
pbX0=(0 1 2 3 4)

#muon runs
for thick in ${siThicknesses[@]}; do 
    for x0 in ${pbX0[@]}; do 
	python scripts/submitFastTimeTestBeamProd.py -q 8nh -n 50000 -e 150 --particle mu- -v ${thick}${x0} -o /store/cmst3/group/hgcal/TimingTB_H2_Jul2015/SIM; 
    done 
done

#electron runs
#for thick in ${siThicknesses[@]}; do 
#    for x0 in ${pbX0[@]}; do 
#	python scripts/submitFastTimeTestBeamProd.py -q 8nh -n 50000 -e 50 --particle e- -v ${thick}${x0} -o /store/cmst3/group/hgcal/TimingTB_H2_Jul2015/SIM; 
#   done 
#done

