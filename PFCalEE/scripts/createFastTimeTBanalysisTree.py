import ROOT
import os
import sys
import optparse
from array import array

"""
Simulation output analysis
"""
def analyze(url,windowSize=2.5,fixedWindow=True,fOutName='H4treeReco.root'):

    #get G4 simulation tree from file and 
    fIn=ROOT.TFile.Open(url)
    tree=fIn.Get('HGCSSTree')
    print '....analyzing %s with %d events' % (url,tree.GetEntriesFast())

    #prepare output file with summary tree
    if fOutName==url :
        print 'Error output file name matches input file name, bailing out...'
        return
    fOut=ROOT.TFile.Open(fOutName,'RECREATE')
    fOut.cd()

    #init summary tree
    outT=ROOT.TTree( 'H4treeReco', 'H4treeReco')
    outT.SetDirectory(fOut)
    runNumber = array( 'i', [ 1 ] )
    outT.Branch( 'runNumber', runNumber, 'runNumber/I' )
    spillNumber = array( 'i', [ 1 ] )
    outT.Branch( 'spillNumber', spillNumber, 'spillNumber/I' )
    evtNumber = array( 'i', [ 0 ] )
    outT.Branch( 'evtNumber', evtNumber, 'evtNumber/I' )
    wc_recox = array( 'f', 2*[ 0 ] )
    outT.Branch( 'wc_recox', wc_recox, 'wc_recox/F' )
    wc_recoy = array( 'f', 2*[ 0 ] )
    outT.Branch( 'wc_recoy', wc_recoy, 'wc_recoy/F' )
    maxch = 10
    nch = array( 'i', [ 0 ] )
    outT.Branch( 'maxch', nch, 'maxch/I' )
    group = array( 'i', maxch*[ 1 ] )
    outT.Branch( 'group', group, 'group[maxch]/I' )
    ch = array( 'i', maxch*[ 0 ] )
    outT.Branch( 'ch', ch, 'ch[maxch]/I' )
    charge_integ = array( 'f', maxch*[ 0. ] )
    outT.Branch( 'charge_integ', charge_integ, 'charge_integ[maxch]/F' )
    charge_integ_full = array( 'f', maxch*[ 0. ] )
    outT.Branch( 'charge_integ_full', charge_integ_full, 'charge_integ_full[maxch]/F' )
    charge_pmult = array( 'f', maxch*[ 0. ] )
    outT.Branch( 'charge_pmult',charge_pmult, 'charge_pmult[maxch]/F' )
    dR_avg = array( 'f', maxch*[ 0. ] )
    outT.Branch( 'dR_avg', dR_avg, 'dR_avg[maxch]/F' )
    dR_rms = array( 'f', maxch*[ 0. ] )
    outT.Branch( 'dR_rms', dR_rms, 'dR_rms[maxch]/F' )

    #loop over events
    for iev in xrange(0,tree.GetEntriesFast()):
        tree.GetEntry(iev)

        #incoming beam position
        evHeader=tree.HGCSSEvent

        evtNumber[0]=evHeader.eventNumber()
        x0,y0=evHeader.vtx_x(),evHeader.vtx_y()
        wc_recox[0]=x0
        wc_recoy[0]=y0
        wc_recox[1]= 0 if fixedWindow else x0
        wc_recoy[1]= 0 if fixedWindow else y0

        #integrate total energy in sim hits
        totalEn,totalEnInWindow={},{}
        partMultInWindow={}
        enWgtR,enWgtR2={},{}
        simhits=tree.HGCSSSimHitVec
        for hit in simhits:

            x=hit.get_x()
            dx=x if fixedWindow else x-x0
            
            y=hit.get_y()
            dy=y if fixedWindow else y-y0

            rho=ROOT.TMath.Sqrt(dx**2+dy**2)
            acceptInWindow=False
            if ROOT.TMath.Abs(dx)<windowSize and ROOT.TMath.Abs(dy)<windowSize: 
                acceptInWindow=True
            
            lay=hit.layer()

            #init counters if not yet started
            if not lay in totalEn: 
                totalEn[lay]=0.0
                totalEnInWindow[lay]=0.0
                partMultInWindow[lay]=0
                enWgtR[lay]=0.0
                enWgtR2[lay]=0.0

            en=hit.energy()*1000.
            partMultInWindow[lay]+=hit.numberOfParticles()
            totalEn[lay]+=en
            enWgtR[lay]+=en*(rho)
            enWgtR2[lay]+=en*(rho**2)
            if acceptInWindow:
                totalEnInWindow[lay]+=en

        #add to summary tree
        maxch=0
        for lay in totalEn:
            ch[maxch]=lay
            charge_pmult[maxch]=partMultInWindow[lay]
            charge_integ[maxch]=totalEnInWindow[lay]
            charge_integ_full[maxch]=totalEn[lay]
            if totalEn[lay]>0:
                dR_avg[maxch]=enWgtR[lay]/totalEn[lay]
                dR_rms[maxch]=ROOT.TMath.Sqrt( (enWgtR2[lay]-(dR_avg[lay]**2))/totalEn[lay] )
            else:
                dR_avg[maxch]=0
                dR_rms[maxch]=0
            maxch+=1
        nch[0]=maxch
        outT.Fill()

    fIn.Close()

    #save tree to file
    fOut.cd()
    outT.Write()
    fOut.Close()

"""
Wrapper to be used when run in parallel
"""
def analyzePacked(args):
    url,windowSize,fixedWindow,fOutName = args
    try:
        analyze(url=url,windowSize=windowSize,fixedWindow=fixedWindow,fOutName=fOutName)
    except :
        print 50*'<'
        print "  Problem  (%s) with %s continuing without"%(sys.exc_info()[1],url)
        print 50*'<'
        return False
    return True



"""
steer the script
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--inDir',   dest='inDir',  help='input directory with files', default='/store/cmst3/group/hgcal/TimingTB_H2_Jul2015/SIM/11ed298', type='string')
    parser.add_option('-o', '--outDir',  dest='outDir', help='output directory',           default='FTTBAnalysis',   type='string')
    parser.add_option('-w', '--window',  dest='window', help='Window size (mm)',           default=2.5,  type=float)
    parser.add_option('-n', '--njobs',   dest='njobs',  help='# jobs to run in parallel',  default=0,    type=int)
    (opt, args) = parser.parse_args()

    #load userlib classes
    ROOT.gSystem.Load('libPFCalEEuserlib')

    #create list of tasks
    task_list=[]

    #mount eos if required
    useEOS=True if '/store/cmst3/' in opt.inDir else False
    if useEOS:
        print 'Mounting eos locally'
        os.system('mkdir -p eos')
        os.system('eos.select -b fuse mount eos')
        opt.inDir = 'eos/cms/'+opt.inDir

    for (dirpath, dirnames, filenames) in os.walk(opt.inDir):
        if len(filenames)==0 : continue
        for f in filenames:
            if not 'root' in f: continue
            url=os.path.join(dirpath,f)
            fOutName=url.replace(opt.inDir,'')
            fOutName=fOutName.replace('/','_')
            fOutName=os.path.join(opt.outDir,fOutName)
            if useEOS: url='root://eoscms//'+url
            fixedWindow=True #False if 'mu-' in url else True
            task_list.append( (url,opt.window,fixedWindow,fOutName) )

    #un-mount
    if useEOS:
        print 'Unmounting eos'
        os.system('eos.select -b fuse umount eos')

    #prepare output
    if len(opt.outDir)!=0 and not os.path.exists(opt.outDir) : os.system( 'mkdir -p %s' % opt.outDir )

    #run the analysis jobs
    if opt.njobs == 0:
        for url,windowSize,fixedWindow,fOutName in task_list:
            analyze(url=url,windowSize=windowSize,fixedWindow=fixedWindow,fOutName=fOutName)
    else:
        from multiprocessing import Pool
        pool = Pool(opt.njobs)
        pool.map(analyzePacked, task_list)
    
"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
