#!/usr/bin/env python

import os,sys
import optparse
import commands
import math
import random

random.seed()

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-s', '--short-queue',    dest='squeue'             , help='short batch queue'            , default='1nd')
parser.add_option('-q', '--long-queue' ,    dest='lqueue'             , help='long batch queue'             , default='2nw')
parser.add_option('-t', '--git-tag'     ,    dest='gittag'             , help='git tag version'              , default='V00-00-00')
parser.add_option('-r', '--run'         ,    dest='run'                , help='stat run'                     , default=-1,      type=int)
parser.add_option('-v', '--version'     ,    dest='version'            , help='detector version'             , default=3,      type=int)
parser.add_option('-m', '--model'       ,    dest='model'              , help='detector model'               , default=3,      type=int)
parser.add_option('-a', '--alpha'       ,    dest='alpha'              , help='incidence angle in rad'       , default=0,      type=float)
parser.add_option('-b', '--Bfield'      ,    dest='Bfield'             , help='B field value in Tesla'       , default=0,      type=float)
parser.add_option('-d', '--datatype'    ,    dest='datatype'           , help='data type or particle to shoot', default='e-')
parser.add_option('-f', '--datafile'    ,    dest='datafile'           , help='full path to HepMC input file', default='data/example_MyPythia.dat')
parser.add_option('-n', '--nevts'       ,    dest='nevts'              , help='number of events to generate' , default=1000,    type=int)
parser.add_option('-o', '--out'         ,    dest='out'                , help='output directory'             , default=os.getcwd() )
parser.add_option('-e', '--eos'         ,    dest='eos'                , help='eos path to save root file to EOS',         default='')
parser.add_option('-g', '--gun'         ,    action="store_true",  dest='dogun'              , help='use particle gun.')
parser.add_option('-S', '--no-submit'   ,    action="store_true",  dest='nosubmit'           , help='Do not submit batch job.')
(opt, args) = parser.parse_args()

enlist=[0]
#if opt.dogun : enlist=[10,15,18,20,25,30,35,40,45,50,60,80,100,200,300,400,500]
#if opt.dogun : enlist=[10,20,30,40,50,100]
if opt.dogun : enlist=[10,15,18,20,25,30,35,40,45,50,60,80] #,100,150,200,300,400,500]

granularity='0-20:4,21-63:6'
noise='0-63:0.12'
threshold='0-63:25'

if opt.version<20 :
    granularity='0-20:4,21-30:6'
    noise='0-30:0.12'
    threshold='0-30:25'
elif (opt.version==21) :
    granularity='0-32:6'
    noise='0-32:0.12'
    threshold='0-32:25'
elif (opt.version==22) :
    granularity='0-8:6'
    noise='0-8:0.12'
    threshold='0-8:25'
elif (opt.version==23) :
    granularity='0-53:12'
    noise='0-53:0.12'
    threshold='0-53:25'

suffix='' #highnoise'

for en in enlist :

    nevents=opt.nevts
    
    myqueue=opt.lqueue
    
    bval="BOFF"
    if opt.Bfield>0 : bval="BON" 
    
    outDir='%s/git_%s/version_%d/model_%d/%s/%s'%(opt.out,opt.gittag,opt.version,opt.model,opt.datatype,bval)
    if en>0 : outDir='%s/e_%d'%(outDir,en)
    eosDir='%s/git%s/%s'%(opt.eos,opt.gittag,opt.datatype)
    if opt.alpha>0 : outDir='%s/a_%3.3f/'%(outDir,opt.alpha) 
    if (opt.run>=0) : outDir='%s/run_%d/'%(outDir,opt.run)

    outlog='%s/digitizer%s.log'%(outDir,suffix)
    g4log='digijob%s.log'%(suffix)
    os.system('mkdir -p %s'%outDir)
    
    #wrapper
    scriptFile = open('%s/runDigiJob%s.sh'%(outDir,suffix), 'w')
    scriptFile.write('#!/bin/bash\n')
    scriptFile.write('source %s/../g4env.sh\n'%(os.getcwd()))
    #scriptFile.write('cd %s\n'%(outDir))
    outTag='version%d_model%d_%s'%(opt.version,opt.model,bval)
    if en>0 : outTag='%s_e%d'%(outTag,en)
    if opt.alpha>0 : outTag='%s_alpha%3.3f'%(outTag,opt.alpha) 
    if (opt.run>=0) : outTag='%s_run%d'%(outTag,opt.run)
    scriptFile.write('localdir=`pwd`\n')
    scriptFile.write('%s/bin/digitizer 0 root://eoscms//eos/cms%s/HGcal_%s.root $localdir/ %s %s %s | tee %s\n'%(os.getcwd(),eosDir,outTag,granularity,noise,threshold,outlog))
    scriptFile.write('echo "--Local directory is " $localdir >> %s\n'%(g4log))
    scriptFile.write('ls * >> %s\n'%(g4log))
    if len(opt.eos)>0:
        scriptFile.write('grep "alias eos=" /afs/cern.ch/project/eos/installation/cms/etc/setup.sh | sed "s/alias /export my/" > eosenv.sh\n')
        scriptFile.write('source eosenv.sh\n')
        scriptFile.write('cmsStage -f DigiPFcal.root %s/Digi%s_%s.root\n'%(eosDir,suffix,outTag))
        scriptFile.write('if (( "$?" != "0" )); then\n')
        scriptFile.write('echo " --- Problem with copy of file DigiPFcal.root to EOS. Keeping locally." >> %s\n'%(g4log))
        scriptFile.write('else\n')
        scriptFile.write('eossize=`$myeos ls -l %s/Digi%s_%s.root | awk \'{print $5}\'`\n'%(eosDir,suffix,outTag))
        scriptFile.write('localsize=`ls -l DigiPFcal.root | awk \'{print $5}\'`\n')
        scriptFile.write('if (( "$eossize" != "$localsize" )); then\n')
        scriptFile.write('echo " --- Copy of digi file to eos failed. Localsize = $localsize, eossize = $eossize. Keeping locally..." >> %s\n'%(g4log))
        scriptFile.write('else\n')
        scriptFile.write('echo " --- Size check done: Localsize = $localsize, eossize = $eossize" >> %s\n'%(g4log))
        scriptFile.write('echo " --- File DigiPFcal.root successfully copied to EOS: %s/Digi%s_%s.root" >> %s\n'%(eosDir,suffix,outTag,g4log))
        scriptFile.write('rm DigiPFcal.root\n')
        scriptFile.write('fi\n')
        scriptFile.write('fi\n')
    scriptFile.write('cp * %s/\n'%(outDir))
    scriptFile.write('echo "All done"\n')
    scriptFile.close()
    
    #submit
    os.system('chmod u+rwx %s/runDigiJob%s.sh'%(outDir,suffix))
    if opt.nosubmit : os.system('echo bsub -q %s %s/runDigiJob%s.sh'%(myqueue,outDir,suffix)) 
    else: os.system("bsub -q %s \'%s/runDigiJob%s.sh\'"%(myqueue,outDir,suffix))

