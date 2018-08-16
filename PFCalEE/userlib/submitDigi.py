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
parser.add_option(      '--nPuVtx'      ,    dest='nPuVtx'             , help='pileup scenarios (csv) [%h]',   default='0,140,200')
parser.add_option('-e', '--eos'         ,    dest='eos'                , help='eos path to save root file to EOS',         default='')
parser.add_option('-E', '--eosin'       ,    dest='eosin'              , help='eos path to read input root file from EOS',  default='')
parser.add_option('-g', '--gun'         ,    action="store_true",  dest='dogun'              , help='use particle gun.')
parser.add_option('-S', '--no-submit'   ,    action="store_true",  dest='nosubmit'           , help='Do not submit batch job.')
parser.add_option(      '--enList'      ,    dest='enList'              , help='E_T list to use with gun [%default]', default='5,10,20,30,40,60,80,100,150,200')
parser.add_option('--etamean' , dest='etamean' , help='mean value of eta ring to save' , default=0,  type=float)
parser.add_option('--deta'    , dest='deta'    , help='width of eta ring'              , default=0,  type=float)
parser.add_option('--forceDigiStep'    , dest='forceDigiStep'    , help='force digitization set digis already found in directory [%default]', default=False, action='store_true')
parser.add_option('--inPathPU'    , dest='inPathPU'    , help='input path for PU files (overrides defaults) [%default]'              , default=None,  type='string')

(opt, args) = parser.parse_args()

#1 = hexagons, 2=diamonds, 3=triangles.
#shape=1
#for run in `seq 0 19`; do ./submitDigi.py -s 1nw -q 1nw -g -t testV8 -r $run -v 63 -m 2 -a 2.0 -b 3.8 -d gamma -n 100 -o /afs/cern.ch/work/a/amagnan/public/HGCalTDR/ -e /store/group/dpg_hgcal/comm_hgcal/amagnan/HGCalTDR -E /store/cmst3/group/hgcal/HGCalTDR ; done
#for run in `seq 0 19`; do ./submitDigi.py -s 1nw -q 1nw -g -t testV8 -r $run -v 63 -m 2 -a 1.7 -b 3.8 -d gamma -n 100 -o /afs/cern.ch/work/a/amagnan/public/HGCalTDR/ -e /store/group/dpg_hgcal/comm_hgcal/amagnan/HGCalTDR -E /store/cmst3/group/hgcal/HGCalTDR ; done




nSiLayers=2

enlist=[0]
#if opt.dogun : enlist=[3,5,7,10,20,30,40,50,60,70,80,90,100,125,150,175,200]
if opt.dogun : 
    #enlist=[50]
    #enlist=[5,10,20,30,50,70,100]
    #enlist=[5,10,20,30,40,60,80,100,150,200]
    enlist=[float(x) for x in opt.enList.split(',')]

#if opt.dogun : enlist=[2,5,10,20,40,60,80,100,150,200]#,300,400,500]

#No label: for 100um Si->will use small cell PU production!
#label=''
#for using Si noise values for 300um and large hexagon cells
#label='300u'
#for using Si noise values for 200um and large hexagon cells
label='200u'
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
elif opt.version==63:
    if (label==''):
        INPATHPU="root://eoscms//eos/cms/store/cmst3/group/hgcal/HGCalTDR/gittestV8/MinBiasSmall/"
    else :
        INPATHPU="root://eoscms//eos/cms/store/group/dpg_hgcal/comm_hgcal/amagnan/HGCalTDR/gitV08-01-00/MinBiasLarge/"

if opt.inPathPU:
    INPATHPU='root://eoscms//eos/cms/%s'%opt.inPathPU
    print 'Default value for INPATHPU overriden with',INPATHPU

#nPuVtxlist=[0,140,200]
nPuVtxlist=[int(x) for x in opt.nPuVtx.split(',')]

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
elif (opt.version==30 or opt.version==100 or opt.version==110):
    granularity='0-27:4'
    noise='0-27:0.14'
    threshold='0-27:5'
elif (opt.version==33):
    #granularity='0-27:4,28-39:4,40-51:8'
    granularity='0-27:1,28-39:1,40-51:1'
    noise='0-39:0.,40-51:0.'
    #noise='0-39:0.14,40-51:0.2'
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
elif (opt.version==60):
    granularity='0-27:1'
    noise='0-27:0.12'
    threshold='0-27:5'
elif (opt.version==61):
    granularity='0-39:1'
    noise='0-23:0.12,24-39:0.15'
    threshold='0-39:5'
elif (opt.version==62):
    granularity='0-15:1'
    noise='0-15:0.15'
    threshold='0-15:5'
elif opt.version==63:
    granularity='0-68:1'
    if (label=='200u'):
        noise='0-51:0.13,53-68:0.15'
    elif (label=='300u'):
        noise='0-51:0.07,53-68:0.15'
    else:
        noise='0-51:0.27,53-68:0.15'
    threshold='0-68:5'
elif opt.version==80:
    granularity='0-69:1'
    noise='0-51:0.13,53-69:0.15'    
    threshold='0-69:5'
elif opt.version==81:
    granularity='0-58:1'
    noise='0-45:0.13,46-58:0.15'    
    threshold='0-58:5'
elif opt.version==82:
    granularity='0-59:1'
    noise='0-45:0.13,46-59:0.15'    
    threshold='0-59:5'
elif opt.version==83:
    #indices[4]-1  indices[4]+1
    granularity='0-76:1'
    noise='0-56:0.13,57-76:0.15'    
    threshold='0-76:5'

else:
    granularity='0-51:4'
    noise='0-51:0.15'
    threshold='0-51:5'

missingFiles=[]
for nPuVtx in nPuVtxlist:

    for interCalib in interCalibList:
        if nPuVtx>0 :
            suffix='Pu%d_IC%d'%(nPuVtx,interCalib)
            myqueue=opt.lqueue
        else :
            suffix='IC%d'%(interCalib)
            myqueue=opt.squeue

        if opt.model!=2 : suffix='%s_Si%d'%(suffix,nSiLayers)
        if (opt.etamean>1.3): suffix='%s_eta%3.2f_%3.2f'%(suffix,opt.etamean-opt.deta,opt.etamean+opt.deta)

        for en in enlist :
            
            bval="BOFF"
            if opt.Bfield>0 : bval="BON" 
            
            outDir='%s/git_%s/version_%d/model_%d/%s/%s'%(opt.out,opt.gittag,opt.version,opt.model,opt.datatype,bval)
            outDir='%s/%s'%(outDir,label) 
            if en>0 : outDir='%s/et_%d'%(outDir,en)

            #eosDirIn='%s'%(opt.eosin)
            if opt.alpha>0 : outDir='%s/eta_%3.3f'%(outDir,opt.alpha) 
            if opt.phi!=0.5 : outDir='%s/phi_%3.3fpi'%(outDir,opt.phi) 
            if (opt.run>=0) : outDir='%s/run_%d'%(outDir,opt.run)
        
            if len(opt.eos)>0:
                eosDir='/eos/cms%s/git%s/%s'%(opt.eos,opt.gittag,opt.datatype)
                eosDirIn='root://eoscms//eos/cms%s/git%s/%s'%(opt.eosin,opt.gittag,opt.datatype)
            else:
                eosDir='%s/'%(outDir)
                eosDirIn='%s/'%(outDir)

            outlog='%s/digitizer%s.log'%(outDir,suffix)
            g4log='digijob%s.log'%(suffix)
            os.system('mkdir -p %s'%outDir)
            
            #wrapper
            scriptFile = open('%s/runDigiJob%s.sh'%(outDir,suffix), 'w')
            scriptFile.write('#!/bin/bash\n')
            scriptFile.write('source %s/../g4env.sh\n'%(os.getcwd()))
            #scriptFile.write('cd %s\n'%(outDir))
            outTag='_version%d_model%d_%s'%(opt.version,opt.model,bval)
            if en>0 : outTag='%s_et%d'%(outTag,en)
            if opt.alpha>0 : outTag='%s_eta%3.3f'%(outTag,opt.alpha) 
            if opt.phi!=0.5 : outTag='%s_phi%3.3fpi'%(outTag,opt.phi) 
            if (opt.run>=0) : outTag='%s_run%d'%(outTag,opt.run)

            simFile='%s/HGcal_%s.root'%(eosDirIn,outTag)
            if not os.path.isfile(simFile.replace('root://eoscms/','')):
                missingFiles.append(simFile)
                continue
            digiFile='%s/Digi%s_%s%s.root'%(eosDir,suffix,label,outTag)
            if os.path.isfile(digiFile) and not opt.forceDigiStep:
                print 'Digi',digiFile,'is already produced but forceProd=False, passing'
                continue
            #continue


            scriptFile.write('localdir=`pwd`\n')
            if (opt.etamean>1.3):
                 scriptFile.write('%s/bin/digitizer %d %s $localdir/ %s %s %s %d %d %d %s %3.2f %3.2f | tee %s\n'%(os.getcwd(),opt.nevts,simFile,granularity,noise,threshold,interCalib,nSiLayers,nPuVtx,INPATHPU,opt.etamean,opt.deta,outlog))
            else:
                scriptFile.write('%s/bin/digitizer %d %s $localdir/ %s %s %s %d %d %d %s | tee %s\n'%(os.getcwd(),opt.nevts,simFile,granularity,noise,threshold,interCalib,nSiLayers,nPuVtx,INPATHPU,outlog))
            scriptFile.write('echo "--Local directory is " $localdir >> %s\n'%(g4log))
            scriptFile.write('ls * >> %s\n'%(g4log))
            if len(opt.eos)>0:
                #scriptFile.write('grep "alias eos=" /afs/cern.ch/project/eos/installation/cms/etc/setup.sh | sed "s/alias /export my/" > eosenv.sh\n')
                #scriptFile.write('source eosenv.sh\n')
                scriptFile.write('eos mkdir -p %s\n'%eosDir)
                scriptFile.write('eos cp $localdir/DigiPFcal.root %s\n'%digiFile)
                scriptFile.write('if (( "$?" != "0" )); then\n')
                scriptFile.write('echo " --- Problem with copy of file DigiPFcal.root to EOS. Keeping locally." >> %s\n'%(g4log))
                scriptFile.write('else\n')
                scriptFile.write('eossize=`eos ls -l %s | awk \'{print $5}\'`\n'%digiFile)
                scriptFile.write('localsize=`ls -l DigiPFcal.root | awk \'{print $5}\'`\n')
                scriptFile.write('if [ $eossize != $localsize ]; then\n')
                scriptFile.write('echo " --- Copy of digi file to eos failed. Localsize = $localsize, eossize = $eossize. Keeping locally..." >> %s\n'%(g4log))
                scriptFile.write('else\n')
                scriptFile.write('echo " --- Size check done: Localsize = $localsize, eossize = $eossize" >> %s\n'%(g4log))
                scriptFile.write('echo " --- File DigiPFcal.root successfully copied to EOS: %s" >> %s\n'%(digiFile,g4log))
                scriptFile.write('rm DigiPFcal.root\n')
                scriptFile.write('fi\n')
                scriptFile.write('fi\n')
            else:
                scriptFile.write('mv DigiPFcal.root Digi%s_%s%s.root\n'%(suffix,label,outTag))

            scriptFile.write('echo "--deleting core files: too heavy!!"\n')
            scriptFile.write('rm core.*\n')
            scriptFile.write('cp * %s/\n'%(outDir))
            scriptFile.write('echo "All done"\n')
            scriptFile.close()
            
            #submit
            os.system('chmod u+rwx %s/runDigiJob%s.sh'%(outDir,suffix))
            if opt.nosubmit : os.system('echo bsub -q %s %s/runDigiJob%s.sh'%(myqueue,outDir,suffix)) 
            else: os.system("bsub -q %s \'%s/runDigiJob%s.sh\'"%(myqueue,outDir,suffix))

missingFiles=set(missingFiles)
if len(missingFiles)>0:
    print 'The following files are missing in action, jobs were not submitted'
    for f in missingFiles:
        print f

