#!/bin/bash

export USERBASE=`pwd`
export G4Build=${USERBASE}/g4build

source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.22.06/x86_64-centos7-gcc48-opt/bin/thisroot.sh 

export BASEINSTALL=/cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-gcc8-opt/
source ${BASEINSTALL}/setup.sh

export FASTJET_INSTALL=/cvmfs/sft.cern.ch/lcg/releases/fastjet/3.3.4-0d9d5/x86_64-centos7-gcc8-opt

export G4BASE=/cvmfs/geant4.cern.ch/geant4/10.7/x86_64-centos7-gcc8-optdeb

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${BASEINSTALL}/lib64:${USERBASE}/userlib/lib:${USERBASE}/analysis/lib

#slc6 setup
#ARCH=x86_64-slc6-gcc46-opt
#source /afs/cern.ch/sw/lcg/contrib/gcc/4.6/${ARCH}/setup.sh 
#export QTHOME=/afs/cern.ch/sw/lcg/external/qt/4.8.4/${ARCH}/
#export G4BASE=/afs/cern.ch/sw/lcg/external/geant4
#export DAWNHOME=/afs/cern.ch/sw/lcg/external/dawn/3_88a/x86_64-slc5-gcc43-opt/
#export G4DAWNFILE_DEST_DIR=${USERBASE}/DawnFiles/
#export HEPMC_DIR=/afs/cern.ch/sw/lcg/external/HepMC/2.06.08/${ARCH}/
#export FASTJET_INSTALL=/afs/cern.ch/sw/lcg/external/fastjet/3.0.3/${ARCH}/

#export BASEINSTALL=/cvmfs/sft.cern.ch/lcg/views/LCG_97/x86_64-centos7-gcc8-opt/

#source ${BASEINSTALL}/setup.sh

cd ${G4BASE}/share/Geant4-10.7.0/geant4make/
source geant4make.sh
cd - &> /dev/null

export G4DIR=${G4BASE}/lib64
mkdir -p $USERBASE/g4build


#export HepMC_DIR=/cvmfs/sft.cern.ch/lcg/releases/HepMC/2.06.10-1a364/x86_64-centos7-gcc8-opt/
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${BASEINSTALL}/lib64:${HepMC_DIR}/lib:$USERBASE/userlib/lib:$USERBASE/analysis/lib



#slc6 setup
#source /afs/cern.ch/sw/lcg/contrib/gcc/4.6/x86_64-slc6/setup.sh
#cd /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.18/${ARCH}/root/
#cd /cvmfs/sft.cern.ch/lcg/releases/ROOT/v6-20-00-patches-5b35b/${ARCH}/
#source bin/thisroot.sh
#cd - &> /dev/null
#export PATH=$DAWNHOME/bin:$PATH:$FASTJET_INSTALL/bin
#export PATH=$PATH:${BASEINSTALL}/bin:${G4Build}
