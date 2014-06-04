#!/bin/sh

INPATH=root://eoscms//eos/cms/store/user/amagnan/HGCalEEGeant4/gitV00-00-02/

for v in 8
do
    for p in e-
    do
	MYDIR=PLOTS/version$v/$p/
	mkdir -p $MYDIR
	for r in 3
	do
	    #for e in 10 15 20 25 30 40
	    #for e in 50 60 80 100 150 200 300 400 500
	    for e in 10 15 20 25 30 40 50 60 80 100 150 200 300 500
	    do
		echo "Processing v${v}_p${p}_e${e}_run${r}"
		./bin/validation 0 $INPATH/$p/ version${v}_model3_BOFF_e${e}_run${r}.root $MYDIR/validation_e${e}_run${r}_300um.root 3 >& v${v}_p${p}_e${e}_run${r}_300um.log &
	    done
	done
    done
done
