#!/usr/bin/env python

import os,sys
import optparse
import commands
import math

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-q', '--queue'      ,    dest='queue'              , help='batch queue'                  , default='1nh')
parser.add_option('-g', '--granularity',    dest='granularity'        , help='granularities \"layer_i-layer_j:factor,layer:factor,...\"'             , default='0-29:4')
parser.add_option('-n', '--noise'      ,    dest='noise'              , help='noise (in Mips) \"layer_i-layer_j:factor,layer:factor,...\"'           , default='0-29:0.1')
parser.add_option('-t', '--thresh'     ,    dest='threshold'          , help='threshold (in ADC) \"layer_i-layer_j:factor,layer:factor,...\"'            , default='0-29:25')
parser.add_option('-N', '--Nevts'      ,    dest='Nevts'              , help='number of events to process' , default=0,    type=int)
parser.add_option('-m', '--miptoadc'   ,    dest='MipToADC'           , help='adc counts per MIP' , default=50,    type=int)
parser.add_option('-r', '--randomseed' ,    dest='seed'               , help='random seed' , default=0,    type=int)
parser.add_option('-d', '--debug'      ,    dest='debug'              , help='debug output' , default=0,    type=int)
parser.add_option('-s', '--scenario'   ,    dest='scenario'           , help='Integer describing scenario.' , default=0,    type=int)
parser.add_option('-o', '--out'        ,    dest='out'                , help='output directory'             , default=os.getcwd() )
parser.add_option('-e', '--eos'        ,    dest='eos'                , help='eos path to save root file to EOS', default='root://eoscms//eos/cms/store/user/amagnan/HGCalHEGeant4/run_0')
parser.add_option('-S', '--no-submit'  ,    action="store_true",  dest='nosubmit'           , help='Do not submit batch job.')
(opt, args) = parser.parse_args()

nevents=opt.Nevts

#myg4dir="/afs/cern.ch/work/a/amagnan/public/HGCalEEGeant4/"
#myg4dir="root://eoscms//eos/cms/store/cmst3/group/hgcal/Geant4"
#myg4dir="root://eoscms//eos/cms/store/user/amagnan/HGCalHEGeant4"

for version in [23]:
    #for particle in ["Gamma","GammaPU","PU"]:
    for particle in ["e-","pi-"]: #"SimplePU"]:
        #for eta in [30,35]:
        #for eta in [20,25,30,35]:
        eosDir='%s/%s'%(opt.eos,particle)
        #for en in [5,10,25,40,50,60,80,100,150,200,300,400,500,1000]:
        for en in [10,15,20,25,30,35,40,45,50,60,80]: #,100,200,300,400,500]:
            #for en in [0,1,2,3,4,5,6,7,8,9]:
            #for en in [5,10,25,50,75,100]:
            #for en in [150,200,300,500]:
            
            inDir='%s/HGcal_version%d_e%d.root'%(eosDir,version,en)
            outDir='%s/version_%d/scenario_%d/%s/e_%d/'%(opt.out,version,opt.scenario,particle,en)
            outlog='%s/digitizer.log'%(outDir)
            os.system('mkdir -p %s'%outDir)
            
                #wrapper
            scriptFile = open('%s/runJob.sh'%(outDir), 'w')
            scriptFile.write('#!/bin/bash\n')
            scriptFile.write('source  %s/../g4env.sh\n'%(os.getcwd()))
            scriptFile.write('localdir=`pwd`\n')
            scriptFile.write('echo "--Local directory is " $localdir >> %s\n'%outlog)
            scriptFile.write('%s/bin/digitizer %d %s $localdir %s %s %s %d %d %d | tee %s\n'%(os.getcwd(),opt.Nevts,inDir,opt.granularity,opt.noise,opt.threshold,opt.MipToADC,opt.seed,opt.debug,outlog))
            scriptFile.write('ls * >> %s\n'%outlog)
            if len(opt.eos)>0:
                outTag='version%d_scenario%d_e%d'%(version,opt.scenario,en)
                scriptFile.write('grep "alias eos=" /afs/cern.ch/project/eos/installation/cms/etc/setup.sh | sed "s/alias /export my/" > eosenv.sh\n')
                scriptFile.write('source eosenv.sh\n')
                scriptFile.write('$myeos mkdir -p %s\n'%eosDir)
                scriptFile.write('myfilepath=`ls DigiPFcal*.root`\n')
                scriptFile.write('echo " --- Root file name to copy to EOS is: "$myfilepath>> %s\n'%outlog)
                scriptFile.write('cmsStage -f $myfilepath %s/HGcal_%s.root\n'%(eosDir,outTag))
                scriptFile.write('if (( "$?" != "0" )); then\n')
                scriptFile.write('echo " --- Problem with copy of file $myfilepath to EOS. Keeping locally." >> %s\n'%outlog)
                scriptFile.write('else\n')
                scriptFile.write('echo " --- File $myfilepath successfully copied to EOS: %s/HGcal_%s.root" >> %s\n'%(eosDir,outTag,outlog))
                scriptFile.write('rm $myfilepath\n')
                scriptFile.write('fi\n')
            scriptFile.write('cp * %s/\n'%(outDir))
            scriptFile.write('echo "All done"\n')
            scriptFile.close()
            
        #submit
            os.system('chmod u+rwx %s/runJob.sh'%outDir)
        #
            if opt.nosubmit : os.system('echo bsub -q %s %s/runJob.sh'%(opt.queue,outDir)) 
            else: os.system("bsub -q %s \'%s/runJob.sh\'"%(opt.queue,outDir))
                
