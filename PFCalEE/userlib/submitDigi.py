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
parser.add_option('-p', '--phi'         ,    dest='phi'                , help='incidence phi angle in pi unit' , default=0.5,      type=float)
parser.add_option('-b', '--Bfield'      ,    dest='Bfield'             , help='B field value in Tesla'       , default=0,      type=float)
parser.add_option('-d', '--datatype'    ,    dest='datatype'           , help='data type or particle to shoot', default='e-')
parser.add_option('-f', '--datafile'    ,    dest='datafile'           , help='full path to HepMC input file', default='data/example_MyPythia.dat')
parser.add_option('-n', '--nevts'       ,    dest='nevts'              , help='number of events to generate' , default=1000,    type=int)
parser.add_option('-o', '--out'         ,    dest='out'                , help='output directory'             , default=os.getcwd() )
parser.add_option('-e', '--eos'         ,    dest='eos'                , help='eos path to save root file to EOS',         default='')
parser.add_option('-E', '--eosin'       ,    dest='eosin'              , help='eos path to read input root file from EOS',  default='')
parser.add_option('-g', '--gun'         ,    action="store_true",  dest='dogun'              , help='use particle gun.')
parser.add_option('-S', '--no-submit'   ,    action="store_true",  dest='nosubmit'           , help='Do not submit batch job.')
(opt, args) = parser.parse_args()

enlist=[0]
#if opt.dogun : enlist=[3,5,7,10,20,30,40,50,60,70,80,90,100,125,150,175,200]
if opt.dogun : 
    enlist=[1,2,3,4,5,6,7,8,9,10,12,14,16,18,20]
    #enlist=[3,5,10,30,50,70,100,200]

#if opt.dogun : enlist=[2,5,10,20,40,60,80,100,150,200]#,300,400,500]

label=''
#label='v5_30'
#label='v5_28'
#label='v5_24'
#label='v5_18'

#INPATHPU="root://eoscms//eos/cms/store/user/msun/V12/MinBias/"
INPATHPU="root://eoscms//eos/cms/store/cmst3/group/hgcal/Standalone/V12/MinBias/"

if opt.version==13:
    INPATHPU="root://eoscms//eos/cms/store/cmst3/group/hgcal/Standalone/V13/MinBias/"
elif opt.version==25:
    INPATHPU="root://eoscms//eos/cms/store/cmst3/group/hgcal/Standalone/V25/MinBias/"
elif opt.version==33:
    INPATHPU="root://eoscms//eos/cms/store/cmst3/group/hgcal/Standalone/V33/MinBias/pile/gitV00-03-07/e-/"

#nPuVtxlist=[0,140,200]
nPuVtxlist=[0]

#in %
interCalibList=[3];#0,1,2,3,4,5,10,15,20,50]

granularity='0-29:4,30-65:4'
noise='0-65:0.15'
threshold='0-65:5'

if (opt.version==8) :
    granularity='0-20:4,21-30:6'
    noise='0-30:0.14'
    threshold='0-30:2'
elif opt.version<20 :
    granularity='0-19:4,20-29:4'
    noise='0-29:0.14'
    threshold='0-29:5'
elif (opt.version==21 or opt.version==24):
    granularity='0-23:6,24-33:8'
    noise='0-33:0.14'
    threshold='0-33:2'
elif opt.version==22:
    granularity='0-9:8'
    noise='0-9:0.14'
    threshold='0-9:2'
elif opt.version==23:
    granularity='0-53:12'
    noise='0-53:0.14'
    threshold='0-53:2'
elif (opt.version==25 or opt.version==26):
    granularity='0-29:4,30-41:4,42-53:8'
    noise='0-41:0.14,42-53:0.2'
    threshold='0-53:5'
elif (opt.version==30):
    granularity='0-27:4'
    noise='0-27:0.14'
    threshold='0-27:5'
elif (opt.version==33):
    granularity='0-27:4,28-39:4,40-51:8'
    noise='0-39:0.14,40-51:0.2'
    threshold='0-51:5'
elif (opt.version==27 or opt.version==31):
    granularity='0-11:4,12-23:8'
    noise='0-11:0.14,12-23:0.2'
    threshold='0-23:5'
elif (opt.version==28 or opt.version==32):
    granularity='0-11:8'
    noise='0-11:0.2'
    threshold='0-11:5'
elif (opt.version==34):
    granularity='0-23:4'
    noise='0-23:0.14'
    threshold='0-23:5'
elif (opt.version==36):
    granularity='0-23:4,24-34:4,35-46:8'
    noise='0-34:0.14,35-46:0.2'
    threshold='0-46:5'
elif (opt.version==38):
    granularity='0-10:4,11-22:8'
    noise='0-10:0.14,11-22:0.2'
    threshold='0-22:5'
elif (opt.version==35):
    granularity='0-17:4'
    noise='0-17:0.14'
    threshold='0-17:5'
elif (opt.version==37):
    granularity='0-17:4,18-26:4,27-38:8'
    noise='0-26:0.14,27-38:0.2'
    threshold='0-38:5'
elif (opt.version==39):
    granularity='0-8:4,9-20:8'
    noise='0-8:0.14,9-20:0.2'
    threshold='0-20:5'
else:
    granularity='0-51:4'
    noise='0-51:0.15'
    threshold='0-51:5'

for nPuVtx in nPuVtxlist:

    for interCalib in interCalibList:
        if nPuVtx>0 :
            suffix='Pu%d_IC%d'%(nPuVtx,interCalib)
            myqueue=opt.lqueue
        else :
            suffix='IC%d'%interCalib
            myqueue=opt.squeue
            
        for en in enlist :
            
            bval="BOFF"
            if opt.Bfield>0 : bval="BON" 
            
            outDir='%s/git_%s/version_%d/model_%d/%s/%s'%(opt.out,opt.gittag,opt.version,opt.model,opt.datatype,bval)
            outDir='%s/%s'%(outDir,label) 
            if en>0 : outDir='%s/et_%d'%(outDir,en)
            eosDir='%s/git%s/%s'%(opt.eos,opt.gittag,opt.datatype)
            eosDirIn='%s/git%s/%s'%(opt.eosin,opt.gittag,opt.datatype)
            #eosDirIn='%s'%(opt.eosin)
            if opt.alpha>0 : outDir='%s/eta_%3.3f/'%(outDir,opt.alpha) 
            if opt.phi!=0.5 : outDir='%s/phi_%3.3fpi/'%(outDir,opt.phi) 
            if (opt.run>=0) : outDir='%s/run_%d/'%(outDir,opt.run)
        
            outlog='%s/digitizer%s.log'%(outDir,suffix)
            g4log='digijob%s.log'%(suffix)
            os.system('mkdir -p %s'%outDir)
            
            #wrapper
            scriptFile = open('%s/runDigiJob%s.sh'%(outDir,suffix), 'w')
            scriptFile.write('#!/bin/bash\n')
            scriptFile.write('source %s/../g4env.sh\n'%(os.getcwd()))
            #scriptFile.write('cd %s\n'%(outDir))
            outTag='%s_version%d_model%d_%s'%(label,opt.version,opt.model,bval)
            if en>0 : outTag='%s_et%d'%(outTag,en)
            if opt.alpha>0 : outTag='%s_eta%3.3f'%(outTag,opt.alpha) 
            if opt.phi!=0.5 : outTag='%s_phi%3.3fpi'%(outTag,opt.phi) 
            if (opt.run>=0) : outTag='%s_run%d'%(outTag,opt.run)
            scriptFile.write('localdir=`pwd`\n')
            scriptFile.write('%s/bin/digitizer %d root://eoscms//eos/cms%s/HGcal_%s.root $localdir/ %s %s %s %d %d %s | tee %s\n'%(os.getcwd(),opt.nevts,eosDirIn,outTag,granularity,noise,threshold,interCalib,nPuVtx,INPATHPU,outlog))
            scriptFile.write('echo "--Local directory is " $localdir >> %s\n'%(g4log))
            scriptFile.write('ls * >> %s\n'%(g4log))
            if len(opt.eos)>0:
                scriptFile.write('grep "alias eos=" /afs/cern.ch/project/eos/installation/cms/etc/setup.sh | sed "s/alias /export my/" > eosenv.sh\n')
                scriptFile.write('source eosenv.sh\n')
                scriptFile.write('$myeos mkdir -p %s\n'%eosDir)
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
            scriptFile.write('echo "--deleting core files: too heavy!!"\n')
            scriptFile.write('rm core.*\n')
            scriptFile.write('cp * %s/\n'%(outDir))
            scriptFile.write('echo "All done"\n')
            scriptFile.close()
            
            #submit
            os.system('chmod u+rwx %s/runDigiJob%s.sh'%(outDir,suffix))
            if opt.nosubmit : os.system('echo bsub -q %s %s/runDigiJob%s.sh'%(myqueue,outDir,suffix)) 
            else: os.system("bsub -q %s \'%s/runDigiJob%s.sh\'"%(myqueue,outDir,suffix))



