export USERBASE=`pwd`
#slc6 setup
#ARCH=x86_64-slc6-gcc46-opt
#source /afs/cern.ch/sw/lcg/contrib/gcc/4.6/${ARCH}/setup.sh 
#export QTHOME=/afs/cern.ch/sw/lcg/external/qt/4.8.4/${ARCH}/
#export G4BASE=/afs/cern.ch/sw/lcg/external/geant4
#export DAWNHOME=/afs/cern.ch/sw/lcg/external/dawn/3_88a/x86_64-slc5-gcc43-opt/
#export G4DAWNFILE_DEST_DIR=${USERBASE}/DawnFiles/
#export HEPMC_DIR=/afs/cern.ch/sw/lcg/external/HepMC/2.06.08/${ARCH}/
#export FASTJET_INSTALL=/afs/cern.ch/sw/lcg/external/fastjet/3.0.3/${ARCH}/

ARCH=x86_64-centos7-gcc62-opt
source /cvmfs/sft.cern.ch/lcg/contrib/gcc/6.2/${ARCH}/setup.sh 

export QTHOME=/cvmfs/sft.cern.ch/lcg/releases/qt/4.8.7-0b84e/${ARCH}/
export G4BASE=/cvmfs/sft.cern.ch/lcg/releases/Geant4/
export HEPMC_DIR=/cvmfs/sft.cern.ch/lcg/releases/HepMC/2.06.09-0a23a/${ARCH}

export FASTJET_INSTALL=/cvmfs/sft.cern.ch/lcg/releases/fastjet/3.3.2-c962b/${ARCH}/


cd $G4BASE/10.06-f2040/${ARCH}/share/Geant4-10.6.0/geant4make/
source geant4make.sh
cd - &> /dev/null
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$XERCESCROOT/lib:$HEPMC_DIR/lib:$USERBASE/userlib/lib:$USERBASE/analysis/lib:$FASTJET_INSTALL/lib
#slc6 setup
#source /afs/cern.ch/sw/lcg/contrib/gcc/4.6/x86_64-slc6/setup.sh
#cd /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.18/${ARCH}/root/
cd /cvmfs/sft.cern.ch/lcg/releases/ROOT/v6-20-00-patches-5b35b/${ARCH}/
source bin/thisroot.sh
cd - &> /dev/null
#export PATH=$DAWNHOME/bin:$PATH:$FASTJET_INSTALL/bin
export PATH=$PATH:$FASTJET_INSTALL/bin
