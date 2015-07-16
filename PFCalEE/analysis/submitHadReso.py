#!/usr/bin/env python

import os,sys
import optparse
import commands
import math
import random

random.seed()

#drop=""
#drop='1,3,5,7,10,12,15,18,20,23,25,27'
#drop='1,10,15,25,27,40'
#drop='1,3,5,7,10,12,15,18,20,23,25,27,36,38,40'
#drop='1,3,5,7,10,12,15,18,20,23,25,27,37,39,41'
#drop='1,3,5,7,10,12,15,18,20,23,25,27,32,36,40'
suffix='tp2411'

#enlist=[3,5,10,30,50,70,100,200]
enlist=[70,100]
for et in enlist :
    #os.system('nohup ./bin/hadronResolution -c scripts/HadronConfigAM.cfg --dropLayers=%s --outFileEM=EMcalibration_%s --genEnergyEM=%d --doFHEMcalib=true --doBHEMcalib=true --doFHBHcalib=false --eeSimFiles=root://eoscms//eos/cms//store/cmst3/group/hgcal/HGCalDescop/gitV04-02-02/gamma/HGcal_v5_30_version12_model2_BOFF_et%d_eta2.500 --fhSimFiles=root://eoscms//eos/cms//store/cmst3/group/hgcal/HGCalDescop/gitV04-02-02/gamma/HGcal_v5_30_version27_model2_BOFF_et%d_eta2.500 --bhSimFiles=root://eoscms//eos/cms//store/cmst3/group/hgcal/HGCalDescop/gitV04-02-02/gamma/HGcal_v5_30_version28_model2_BOFF_et%d_eta2.500 --eeRecoFiles=root://eoscms//eos/cms//store/cmst3/group/hgcal/HGCalDescop/gitV04-02-02/gamma/DigiIC3_v5_30_version12_model2_BOFF_et%d_eta2.500 --fhRecoFiles=root://eoscms//eos/cms//store/cmst3/group/hgcal/HGCalDescop/gitV04-02-02/gamma/DigiIC3_v5_30_version27_model2_BOFF_et%d_eta2.500 --bhRecoFiles=root://eoscms//eos/cms//store/cmst3/group/hgcal/HGCalDescop/gitV04-02-02/gamma/DigiIC3_v5_30_version28_model2_BOFF_et%d_eta2.500 >& hadReso_%s_e%d.log & '%(drop,suffix,et,et,et,et,et,et,et,suffix,et))
    #os.system('nohup ./bin/hadronResolution -c scripts/HadronConfigAM.cfg --dropLayers=%s --genEnergy=%d --doFHEMcalib=false --doBHEMcalib=false --doFHBHcalib=false --outFilePi=pion_gc3-50_reso_%s >& hadResoPi_%s_e%d.log & '%(drop,et,suffix,suffix,et))
    #os.system('nohup ./bin/hadronResolution -c scripts/HadronConfigAM.cfg --nRunsEM=4 --outFileEM=EMcalibration_%s --genEnergyEM=%d --doFHEMcalib=true --doBHEMcalib=true --doFHBHcalib=false --outPath=/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalDescop/gitV05-02-04/version33/pi-/ --eeSimFiles=root://eoscms//eos/cms//store/cmst3/group/hgcal/HGCalDescop/gitV05-02-04/gamma/HGcal__version35_model2_BOFF_et%d_eta2.500 --fhSimFiles=root://eoscms//eos/cms//store/cmst3/group/hgcal/HGCalDescop/gitV05-02-04/gamma/HGcal__version39_model2_BOFF_et%d_eta2.500 --bhSimFiles=root://eoscms//eos/cms//store/cmst3/group/hgcal/HGCalDescop/gitV05-02-04/gamma/HGcal__version32_model2_BOFF_et%d_eta2.500 --eeRecoFiles=root://eoscms//eos/cms//store/cmst3/group/hgcal/HGCalDescop/gitV05-02-04/gamma/DigiIC3__version35_model2_BOFF_et%d_eta2.500 --fhRecoFiles=root://eoscms//eos/cms//store/cmst3/group/hgcal/HGCalDescop/gitV05-02-04/gamma/DigiIC3__version39_model2_BOFF_et%d_eta2.500 --bhRecoFiles=root://eoscms//eos/cms//store/cmst3/group/hgcal/HGCalDescop/gitV05-02-04/gamma/DigiIC3__version32_model2_BOFF_et%d_eta2.500 >& hadReso_%s_e%d.log & '%(suffix,et,et,et,et,et,et,et,suffix,et))
    os.system('nohup ./bin/hadronResolution -c scripts/HadronConfigAM.cfg --nRuns=10 --filePath=root://eoscms//eos/cms/store/cmst3/group/hgcal/HGCalDescop/gitV05-02-04/pi-/ --simFileName=HGcal__version36_model2_BOFF_et3_eta2.500_run0.root --recoFileName=DigiIC3__version36_model2_BOFF_et3_eta2.500_run0.root --outPath=/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalDescop/gitV05-02-04/version33/pi-/ --genEnergy=%d --doFHEMcalib=false --doBHEMcalib=false --doFHBHcalib=false --outFilePi=pion_gc3-50_reso_%s >& hadResoPi_%s_e%d.log & '%(et,suffix,suffix,et))
