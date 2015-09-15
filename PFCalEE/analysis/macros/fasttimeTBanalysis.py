import ROOT


def analyze(url,windowSize=1.5,fixedWindow=True):

    #get tree from file
    fIn=ROOT.TFile.Open(url)
    tree=fIn.Get('HGCSSTree')
    print '....analyzing %s with %d events' % (url,tree.GetEntriesFast())

    histos={}

    #loop over events
    for iev in xrange(0,tree.GetEntriesFast()):
        tree.GetEntry(iev)
    
        #incoming beam position
        evHeader=tree.HGCSSEvent
        x0,y0=evHeader.vtx_x(),evHeader.vtx_y()

        #integrate total energy in sim hits
        totalEn,totalEnInWindow={},{}
        enWgtR,enWgtR2={},{}
        simhits=tree.HGCSSSimHitVec
        for hit in simhits:

            x=10.*hit.position().x()
            dx=x if fixedWindow else x-x0
            
            y=10.*hit.position().y()
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
                enWgtR[lay]=0.0
                enWgtR2[lay]=0.0

            en=hit.energy()*1000.
            totalEn[lay]+=en
            enWgtR[lay]+=en*(rho)
            enWgtR2[lay]+=en*(rho**2)
            if acceptInWindow:
                totalEnInWindow[lay]+=en

        #add to histograms
        for lay in totalEn:
            key='_%d' % lay
            if not 'totalEn'+key in histos:
                histos['totalEn'+key]=ROOT.TH1F('totalEn'+key,';Energy [keV];Events',500,0,1000)
                histos['totalEn'+key].SetDirectory(0)
                histos['totalEn'+key].Sumw2()
            
                histos['totalEnInWindow'+key]=ROOT.TH1F('totalEnInWindow'+key,';Energy [keV];Events',500,0,1000)
                histos['totalEnInWindow'+key].SetDirectory(0)
                histos['totalEnInWindow'+key].Sumw2()

                histos['avgdR'+key]=ROOT.TH1F('avgdR'+key,';<#rho> [cm];Events',500,0,3.0)
                histos['avgdR'+key].SetDirectory(0)
                histos['avgdR'+key].Sumw2()

                histos['rmsdR'+key]=ROOT.TH1F('rmsdR'+key,';RMS(#rho) [cm];Events',500,0,3.0)
                histos['rmsdR'+key].SetDirectory(0)
                histos['rmsdR'+key].Sumw2()

            if totalEn[lay]==0 : continue
            histos['totalEn'+key].Fill(totalEn[lay])
            histos['totalEnInWindow'+key].Fill(totalEnInWindow[lay])
            avgdR=enWgtR[lay]/totalEn[lay]
            histos['avgdR'+key].Fill(avgdR)
            rmsdR=ROOT.TMath.Sqrt( (enWgtR2[lay]-(avgdR**2))/totalEn[lay] )
            histos['rmsdR'+key].Fill(rmsdR)

    fIn.Close()
    return histos


#load userlib classes
ROOT.gSystem.Load('libPFCalEEuserlib')
histos=analyze('../../da7e420/version_3203/model_5/TBSim_e-_50.root')
#histos=analyze('../../da7e420/version_3203/model_5/TBSim_mu-_150.root',1.5,False)
fOut=ROOT.TFile.Open('TBspectra.root','RECREATE')
for h in histos:
    histos[h].SetDirectory(fOut)
    histos[h].Write()
fOut.Close()
        
