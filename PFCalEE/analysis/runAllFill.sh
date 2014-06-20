#!/bin/sh

INPATH=root://eoscms//eos/cms/store/user/amagnan/HGCalEEGeant4/gitV00-01-01/

for v in 12
do
    for p in e-
    do
	MYDIR=PLOTS/gitV00-01-01/version$v/$p/
	mkdir -p $MYDIR
	for r in 0
	do
	    #for e in 10 15 18 20 25 30 35 40 45 50 60 80
	    #for e in 50 60 80 100 150 200 300 400 500
	    for e in 5 10 15 20 25 30 40 50 60 80 100 200 300 500
	    do
		echo "Processing v${v}_p${p}_e${e}_run${r}"
		./bin/validation 0 $INPATH/$p/ HGcal_version${v}_model1_BOFF_e${e}.root Digi_version${v}_model1_BOFF_e${e}.root $MYDIR/validation_e${e}.root 4 >& v${v}_p${p}_e${e}.log &
	    done
	done
    done
done
