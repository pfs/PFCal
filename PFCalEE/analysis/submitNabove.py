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
parser.add_option('-b', '--Bfield'      ,    dest='Bfield'             , help='B field value in Tesla'       , default=0,      type=float)
parser.add_option('-d', '--datatype'    ,    dest='datatype'           , help='data type or particle to shoot', default='e-')
parser.add_option('-n', '--nevts'       ,    dest='nevts'              , help='number of events to process' , default=0,    type=int)
parser.add_option('-o', '--out'         ,    dest='out'                , help='output directory'             , default=os.getcwd() )
parser.add_option('-E', '--eosin'       ,    dest='eosin'              , help='eos path to read input root file from EOS',  default='')
parser.add_option('-g', '--gun'         ,    action="store_true",  dest='dogun'              , help='use particle gun.')
parser.add_option('-S', '--no-submit'   ,    action="store_true",  dest='nosubmit'           , help='Do not submit batch job.')
(opt, args) = parser.parse_args()

label='_'

workdir='/afs/cern.ch/work/a/amagnan/PFCalEEAna/'

enlist=[0]
if opt.dogun : 
    #enlist=[3,5,7,10,20,30,40,50,60,70,80,90,100,125,150,175,200]
    #enlist=[3,5,10,30,50,70,100,200]
    #enlist=[1,2,3,4,5,6,7,8,9,10,12,14,16,18,20]
    enlist=[20,50,70]
    #enlist=[3,5,7,10,30,50,70,90,125,150,175,200]

#alphaset=[1.750,2.000,2.250,2.500,2.750]
alphaset=[1.700,1.900,2.100,2.300,2.500,2.700]
#alphaset=[1.750,2.000]
#alphaset=[2.500]

interCalibList=[3] #0,1,2,3,4,5,10,15,20,50]

for alpha in alphaset :
    for et in enlist :
        
        nevents=opt.nevts
        myqueue=opt.lqueue
        bval="BOFF"
        if opt.Bfield>0 : bval="BON" 
        
        #too many files produced: make local output then copy back to afs
        #outDir='%s/%s/git%s/version%d/%s/200um/eta%s_et%s_pu%s'%(os.getcwd(),opt.out,opt.gittag,opt.version,opt.datatype,eta,et,nPuVtx)
        outDir='%s/git%s/%s/v%d_et%s_eta%s/'%(opt.out,opt.gittag,opt.datatype,opt.version,et,alpha)
        eosDirIn='%s/git%s/%s'%(opt.eosin,opt.gittag,opt.datatype)
        
        outlog='nabove.log'
        g4log='nabovejob.log'
        os.system('mkdir -p %s/%s'%(workdir,outDir))
        #clean up old batch outputs
        os.system('rm -f %s/%s/*.*.out'%(workdir,outDir))
        #wrapper
        scriptFile = open('%s/%s/runJob.sh'%(workdir,outDir), 'w')
        scriptFile.write('#!/bin/bash\n')
        scriptFile.write('source %s/../g4env.sh\n'%(os.getcwd()))
        scriptFile.write('localdir=`pwd`\n')
        scriptFile.write('cp -r %s/data .\n'%os.getcwd())
        scriptFile.write('cp -r %s/scripts .\n'%os.getcwd())
        scriptFile.write('mkdir -p %s\n'%outDir)
        scriptFile.write('%s/bin/plotNabove100fC 26 root://eoscms//eos/cms%s/ %d %s %s %d %d | tee %s\n'%(os.getcwd(),eosDirIn,opt.version,et,alpha,opt.nRuns,opt.nevts,outlog))
        scriptFile.write('echo "--Local directory is " $localdir >> %s\n'%(g4log))
        scriptFile.write('ls * >> %s\n'%(g4log))
        scriptFile.write('echo "--deleting core files: too heavy!!" >> %s\n'%(g4log))
        scriptFile.write('rm -f core.* >> %s\n'%(g4log))
        scriptFile.write('cp %s/* %s/%s/\n'%(outDir,workdir,outDir))
        scriptFile.write('cp * %s/%s/\n'%(workdir,outDir))
        scriptFile.write('echo "All done"\n')
        scriptFile.close()
                
        #submit
        os.system('chmod u+rwx %s/%s/runJob.sh'%(workdir,outDir))
        if opt.nosubmit : os.system('echo bsub -q %s %s/%s/runJob.sh'%(myqueue,workdir,outDir)) 
        else: os.system("bsub -q %s \'%s/%s/runJob.sh\'"%(myqueue,workdir,outDir))

