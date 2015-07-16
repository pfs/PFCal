#!/usr/bin/env python

import os,sys
import optparse
import commands
import math
import random

random.seed()

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-s', '--short-queue',    dest='squeue'             , help='short batch queue'            , default='1nh')
parser.add_option('-q', '--long-queue' ,    dest='lqueue'             , help='long batch queue'             , default='1nd')
parser.add_option('-t', '--git-tag'     ,    dest='gittag'             , help='git tag version'              , default='V00-00-00')
parser.add_option('-r', '--run'         ,    dest='run'                , help='stat run'                     , default=-1,      type=int)
parser.add_option('-R', '--nRuns'       ,    dest='nRuns'              , help='number of runs'               , default=0,      type=int)
parser.add_option('-v', '--version'     ,    dest='version'            , help='detector version'             , default=3,      type=int)
parser.add_option('-m', '--model'       ,    dest='model'              , help='detector model'               , default=3,      type=int)
parser.add_option('-a', '--alpha'       ,    dest='alpha'              , help='incidence angle in rad'       , default=0,      type=float)
parser.add_option('-p', '--phi'         ,    dest='phi'                , help='incidence phi angle in pi unit' , default=0.5,      type=float)
parser.add_option('-b', '--Bfield'      ,    dest='Bfield'             , help='B field value in Tesla'       , default=0,      type=float)
parser.add_option('-d', '--datatype'    ,    dest='datatype'           , help='data type or particle to shoot', default='e-')
parser.add_option('-f', '--datafile'    ,    dest='datafile'           , help='full path to HepMC input file', default='data/example_MyPythia.dat')
parser.add_option('-n', '--nevts'       ,    dest='nevts'              , help='number of events to process' , default=0,    type=int)
parser.add_option('-o', '--out'         ,    dest='out'                , help='output directory'             , default=os.getcwd() )
parser.add_option('-e', '--eos'         ,    dest='eos'                , help='eos path to save root file to EOS',         default='')
parser.add_option('-E', '--eosin'       ,    dest='eosin'              , help='eos path to read input root file from EOS',  default='')
parser.add_option('-g', '--gun'         ,    action="store_true",  dest='dogun'              , help='use particle gun.')
parser.add_option('-S', '--no-submit'   ,    action="store_true",  dest='nosubmit'           , help='Do not submit batch job.')
(opt, args) = parser.parse_args()

redofit=1
label='v5_30_'
#label='fixedVtx_'
#label='v5_28_'
#label='v5_24_'
#label='v5_18_'

workdir='/afs/cern.ch/work/a/amagnan/PFCalEEAna/'

enlist=[0]
if opt.dogun : 
    #enlist=[3,5,7,10,20,30,40,50,60,70,80,90,100,125,150,175,200]
    #enlist=[3,5,10,30,50,70,100,200]
    enlist=[70]
    #enlist=[3,5,7,10,30,50,70,90,125,150,175,200]

alphaset=[1.600]
#alphaset=[0.361,0.297,0.244,0.200,0.164,0.134,0.110]
#alphaset=[0.361,0.244,0.164,0.110]
#alphaset=[1.700,1.900,2.100,2.300,2.500,2.700]
#alphaset=[0.297,0.244,0.200,0.134,0.110]
#nPuVtxset=[0,140]
nPuVtxset=[0]
#etaset=[17,19,21,23,25,27,29]
#etaset=[17,21,25,29]
etaset=[16]
#etaset=[17,19,21,23,25,27]

interCalibList=[3] #0,1,2,3,4,5,10,15,20,50]

for nPuVtx in nPuVtxset :
    for interCalib in interCalibList:
        counter=0

        for alpha in alphaset :
            eta=etaset[counter]
            counter=counter+1
            for et in enlist :
                
                nevents=opt.nevts
                myqueue=opt.lqueue
                bval="BOFF"
                if opt.Bfield>0 : bval="BON" 
                
            #too many files produced: make local output then copy back to afs
            #outDir='%s/%s/git%s/version%d/%s/200um/eta%s_et%s_pu%s'%(os.getcwd(),opt.out,opt.gittag,opt.version,opt.datatype,eta,et,nPuVtx)
                outDir='%s/git%s/version%d/%s/%s/eta%s_et%s_pu%s_IC%d'%(opt.out,opt.gittag,opt.version,opt.datatype,label,eta,et,nPuVtx,interCalib)
                if opt.phi!=0.5 : outDir='%s/git%s/version%d/%s/200um/phi_%3.3fpi/eta%s_et%s_pu%s_IC%d'%(opt.out,opt.gittag,opt.version,opt.datatype,opt.phi,eta,et,nPuVtx,interCalib)
                eosDir='%s/git%s/%s'%(opt.eos,opt.gittag,opt.datatype)
                eosDirIn='%s/git%s/%s'%(opt.eosin,opt.gittag,opt.datatype)

                outlog='egreso.log'
                g4log='egresojob.log'
                os.system('mkdir -p %s/%s'%(workdir,outDir))
            #clean up old batch outputs
                os.system('rm -f %s/%s/*.*.out'%(workdir,outDir))
            #clean-up evt-by-evt files
                os.system('rm -f %s/%s/initialPos_*'%(workdir,outDir))
            #wrapper
                scriptFile = open('%s/%s/runEGResoJob.sh'%(workdir,outDir), 'w')
                scriptFile.write('#!/bin/bash\n')
                scriptFile.write('source %s/../g4env.sh\n'%(os.getcwd()))
            #scriptFile.write('cd %s\n'%(outDir))
                outTag='%sversion%d_model%d_%s'%(label,opt.version,opt.model,bval)
                if et>0 : outTag='%s_et%d'%(outTag,et)
                if alpha>0 : outTag='%s_eta%3.3f'%(outTag,alpha) 
                if opt.phi!=0.5 : outTag='%s_phi%3.3fpi'%(outTag,opt.phi) 
                if (opt.run>=0) : outTag='%s_run%d'%(outTag,opt.run)
                scriptFile.write('localdir=`pwd`\n')
                scriptFile.write('cp -r %s/data .\n'%os.getcwd())
                scriptFile.write('cp -r %s/scripts .\n'%os.getcwd())
                scriptFile.write('mkdir -p %s\n'%outDir)
                scriptFile.write('cp %s/%s/*.dat %s/.\n'%(workdir,outDir,outDir))
                if (nPuVtx==0) :
                    if (opt.nRuns==0) :
                        scriptFile.write('%s/bin/egammaResolution -c scripts/DefaultConfig.cfg -n %s -i root://eoscms//eos/cms%s/ --digifilePath=root://eoscms//eos/cms%s/ -s HGcal_%s.root -r Digi_%s.root -o %s.root --redoStep=%s | tee %s\n'%(os.getcwd(),opt.nevts,eosDirIn,eosDir,outTag,outTag,outDir,redofit,outlog))
                        #scriptFile.write('%s/bin/egammaResolution -c scripts/DefaultConfig.cfg -n %s -i root://eoscms//eos/cms%s/ --digifilePath=root://eoscms//eos/cms%s/ -s HGcal_%s.root -r DigiIC%d_%s.root -o %s.root --redoStep=%s | tee %s\n'%(os.getcwd(),opt.nevts,eosDirIn,eosDir,outTag,interCalib,outTag,outDir,redofit,outlog))
                    else:
                        scriptFile.write('%s/bin/egammaResolution -c scripts/DefaultConfig.cfg -n %s --nRuns=%s -i root://eoscms//eos/cms%s/ --digifilePath=root://eoscms//eos/cms%s/ -s HGcal_%s -r DigiIC%d_%s -o %s.root --redoStep=%s | tee %s\n'%(os.getcwd(),opt.nevts,opt.nRuns,eosDirIn,eosDir,outTag,interCalib,outTag,outDir,redofit,outlog))
                else:
                    if (opt.nRuns==0) :
                        scriptFile.write('%s/bin/egammaResolution -c scripts/DefaultConfig.cfg -n %s -i root://eoscms//eos/cms%s/ --digifilePath=root://eoscms//eos/cms%s/ -s HGcal_%s.root -r PuMix%s_%s.root -o %s.root --redoStep=%s | tee %s\n'%(os.getcwd(),opt.nevts,eosDirIn,eosDir,outTag,nPuVtx,outTag,outDir,redofit,outlog)) 
                    else:
                        scriptFile.write('%s/bin/egammaResolution -c scripts/DefaultConfig.cfg -n %s --nRuns=%s -i root://eoscms//eos/cms%s/ --digifilePath=root://eoscms//eos/cms%s/ -s HGcal_%s -r PuMix%s_%s -o %s.root --redoStep=%s | tee %s\n'%(os.getcwd(),opt.nevts,opt.nRuns,eosDirIn,eosDir,outTag,nPuVtx,outTag,outDir,redofit,outlog)) 
                        
                scriptFile.write('echo "--Local directory is " $localdir >> %s\n'%(g4log))
                scriptFile.write('ls * >> %s\n'%(g4log))
                scriptFile.write('echo "--deleting core files: too heavy!!" >> %s\n'%(g4log))
                scriptFile.write('rm -f core.* >> %s\n'%(g4log))
                scriptFile.write('echo "--deleting evt-by-evt files: too many!!" >> %s\n'%(g4log))
                scriptFile.write('rm -f %s/initialPos_* >> %s\n'%(outDir,g4log))
                scriptFile.write('cp %s/* %s/%s/\n'%(outDir,workdir,outDir))
                if (redofit==1):
                    scriptFile.write('cp %s.root %s/%s.root\n'%(outDir,workdir,outDir))
                else:
                    scriptFile.write('cp %s.root %s/%s_nofit.root\n'%(outDir,workdir,outDir))
                scriptFile.write('cp * %s/%s/\n'%(workdir,outDir))
                scriptFile.write('echo "All done"\n')
                scriptFile.close()
                
        #submit
                os.system('chmod u+rwx %s/%s/runEGResoJob.sh'%(workdir,outDir))
                if opt.nosubmit : os.system('echo bsub -q %s %s/%s/runEGResoJob.sh'%(myqueue,workdir,outDir)) 
                else: os.system("bsub -q %s \'%s/%s/runEGResoJob.sh\'"%(myqueue,workdir,outDir))

