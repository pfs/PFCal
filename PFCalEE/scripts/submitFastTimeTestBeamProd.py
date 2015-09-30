#!/usr/bin/env python

import os,sys
import optparse
import commands
import math
import random
import time

seed=random.seed()
githash=commands.getstatusoutput('git log --pretty=format:\'%h\' -n 1')[1]
detectorModel=5

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-q', '--queue'       ,    dest='queue'              , help='batch queue'                  , default='local')
parser.add_option('-t', '--git-tag'     ,    dest='gittag'             , help='git tag version'              , default=githash)
parser.add_option('-v', '--version'     ,    dest='version'            , help='detector version'             , default=3203,    type=int)
parser.add_option(      '--particle'    ,    dest='particle'           , help='particle type'                , default='e-')
parser.add_option('-n', '--nevts'       ,    dest='nevents'            , help='number of events to generate' , default=100,     type=int)
parser.add_option('-e', '--en'          ,    dest='energy'             , help='energy'                       , default=50,      type=float)
parser.add_option('-o', '--out'         ,    dest='out'                , help='output directory'             , default=os.getcwd() )
(opt, args) = parser.parse_args()
    
#prepare output directory
timeTag=time.time()
outDir='%s/%s/version_%d/model_%d/'%(opt.out,opt.gittag,opt.version,detectorModel)
outFile='TBSim_%s_%g'%(opt.particle,opt.energy)
if '/store/' in outDir:
    os.system('cmsMkdir %s' % outDir)
else:
    os.system('mkdir -p %s' % outDir)
os.system('mkdir -p JOBS')

#wrapper to run the job
scriptFile = open('JOBS/runJob_%s.sh'%timeTag, 'w')
scriptFile.write('#!/bin/bash\n')
scriptFile.write('source %s/g4env.sh\n'%(os.getcwd()))
scriptFile.write('cp %s/JOBS/g4steer_%s.mac ./\n'%(os.getcwd(),timeTag))
scriptFile.write('PFCalEE g4steer_%s.mac %d %d 0 | tee %s.log\n'%(timeTag,opt.version,detectorModel,outFile))
scriptFile.write('mv PFcal.root %s.root\n' % outFile)
scriptFile.write('localdir=`pwd`\n')
scriptFile.write('echo "--Local directory is $localdir" >> %s.log\n'%outFile)
scriptFile.write('ls * >> %s.log\n'%outFile)
if '/store/' in outDir:
    scriptFile.write('cmsStage %s.root %s/\n' % (outFile,outDir) )
    scriptFile.write('cmsStage %s.log %s/\n' % (outFile,outDir) )
    scriptFile.write('rm %s.*\n' % outFile)
elif outDir!=os.getcwd():
    scriptFile.write('mv %s.root %s/\n' % (outFile,outDir) )
    scriptFile.write('mv %s.log %s/\n' % (outFile,outDir) )
    scriptFile.write('rm %s.*\n' % outFile)
scriptFile.write('rm core.*\n')
scriptFile.write('rm g4steer_%s.mac\n' % timeTag)
scriptFile.write('echo "All done"\n')
scriptFile.close()
    
#write geant 4 macro
g4Macro = open('JOBS/g4steer_%s.mac' % timeTag, 'w')
g4Macro.write('/control/verbose 0\n')
g4Macro.write('/control/saveHistory\n')
g4Macro.write('/run/verbose 0\n')
g4Macro.write('/event/verbose 0\n')
g4Macro.write('/tracking/verbose 0\n')
g4Macro.write('/N03/det/setField 0 T\n')
g4Macro.write('/N03/det/setModel %d\n'%detectorModel)
g4Macro.write('/random/setSeeds %d %d\n'%( random.uniform(0,100000), random.uniform(0,100000) ) )
g4Macro.write('/generator/select particleGun\n')
g4Macro.write('/gun/particle %s\n' % opt.particle )
g4Macro.write('/gun/energy %f GeV\n' % opt.energy )
g4Macro.write('/gun/direction 0.0 0.0 1.0\n')
g4Macro.write('/run/beamOn %d\n'%(opt.nevents))
g4Macro.close()
os.system('chmod u+rwx JOBS/runJob_%s.sh'%timeTag)
if opt.queue != 'local':
    os.system('bsub -q %s %s/JOBS/runJob_%s.sh' % (opt.queue,os.getcwd(),timeTag))
else: 
    print 'Local queue: running job locally'
    os.system('sh JOBS/runJob_%s.sh' % timeTag)

