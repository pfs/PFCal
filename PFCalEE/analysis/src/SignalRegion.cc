#include "SignalRegion.hh"
#include "HGCSSEvent.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGenParticle.hh"

SignalRegion::SignalRegion(const std::string inputFolder,
			   const unsigned nLayers,
			   const unsigned nevt,
			   const HGCSSGeometryConversion & geomConv,
			   const HGCSSPUenergy & puDensity,
			   const bool applyPuMixFix,
			   const unsigned versionNumber){

  nSR_ = 5;
  nevt_ = nevt;
  inputFolder_ = inputFolder;
  nLayers_ = nLayers;
  geomConv_ = geomConv;
  puDensity_ = puDensity;
  fixForPuMixBug_ = applyPuMixFix;
  firstEvent_ = true;
  nSkipped_ = 0;

  double zpos(0);
  int layerIndex(0);
  
  std::ifstream fzpos;
  std::ostringstream finname;
  finname << "data/zPositions_v" << versionNumber << ".dat";
  fzpos.open(finname.str());
  if (!fzpos.is_open()){
    std::cout << " Cannot open input file " << finname.str() << "! Exiting..." << std::endl;
    exit(1);
  }
  else {
    std::cout << " -- z positions taken from file " << finname.str() << std::endl;
  }
  for(unsigned iL(0);iL<nLayers_;iL++){
    fzpos >> layerIndex >> zpos;
    zPos_.push_back(zpos);
  }
}

SignalRegion::~SignalRegion(){
}

bool SignalRegion::initialiseFitPositions(){

  std::ifstream fxypos;
  std::ostringstream finname;
  finname << inputFolder_ << "/accuratePos.dat";
  fxypos.open(finname.str());
  if (!fxypos.is_open()){
    std::cout << " Cannot open input file " << finname.str() << "! Exiting..." << std::endl;
    return false;
  }

  std::cout << " -- Accurate positions found in file " << finname.str() << std::endl;
  
  //all events did not pass the chi2 fit: fill only those found.
  //keep failed ones to emptyvec so they are ignored afterwards
  accurateFit_.clear();
  FitResult res;
  accurateFit_.resize(nevt_,res);

  unsigned nfound = 0;

  while (!fxypos.eof()){
    unsigned eventIndex = nevt_;
    double xpos(0),ypos(0),xangle(0),yangle(0);
    std::string fitMatrix[4] = {"0","0","0","0"};
    double xpost(0),ypost(0),xanglet(0),yanglet(0);
    fxypos >> eventIndex >> xpos >> fitMatrix[0] >> xangle >> fitMatrix[1] >> ypos >> fitMatrix[2] >> yangle >> fitMatrix[3] >> xpost >> xanglet >> ypost >>yanglet ;
    //testing for nan
    //if ( eventIndex != eventIndex || xpos != xpos || fitMatrix[0]!=fitMatrix[0] || xangle!=xangle || fitMatrix[1]!=fitMatrix[1] || ypos!=ypos || fitMatrix[2]!=fitMatrix[2] || yangle!=yangle || fitMatrix[3]!=fitMatrix[3]){
    //std::cout << " Found nan ! Fix code !" << std::endl;
    //std::cout << eventIndex << " " << xpos << " " << fitMatrix[0] << " " << xangle << " " << fitMatrix[1] << " " << ypos << " " << fitMatrix[2] << " " << yangle << " " << fitMatrix[3]<< std::endl;
    //exit(1);
    //}
    if (eventIndex<nevt_) {
      accurateFit_[eventIndex].pos_x = xpos;
      accurateFit_[eventIndex].pos_y = ypos;
      accurateFit_[eventIndex].tanangle_x = xangle;
      accurateFit_[eventIndex].tanangle_y = yangle;
      accurateFit_[eventIndex].found=true;
      nfound++;
    }
    else break;
  }
  //if not all events found
  if (nfound < nevt_) {
    std::cout << " Warning, file " << finname.str() << " contains only " << nfound 
	      << " events, program running on " << nevt_ 
	      << "." << std::endl;
    if (nfound*1./nevt_ < 0.5) return false;
  }
  std::cout << " -- Now filling signal region histograms..." << std::endl;
  return true;

}

void SignalRegion::initialise(TFile *outputFile,
			      const std::string outputDir){

  //mycalib_ = mycalib;
  outputDir_ = outputDir;
  setOutputFile(outputFile);
  initialiseHistograms();
}

const FitResult & SignalRegion::getAccurateFit(const unsigned ievt) const{
  return accurateFit_[ievt];
}

ROOT::Math::XYZPoint SignalRegion::getAccuratePos(const unsigned ievt, const unsigned iL) const{
  return getAccuratePos(accurateFit_[ievt],iL);
}

ROOT::Math::XYZPoint SignalRegion::getAccuratePos(const FitResult& fit, const unsigned iL) const{
  return ROOT::Math::XYZPoint(fit.pos_x + fit.tanangle_x*zPos_[iL], fit.pos_y + fit.tanangle_y*zPos_[iL], zPos_[iL]);
}

Direction SignalRegion::getAccurateDirection(const unsigned ievt) const{
  return Direction(accurateFit_[ievt].tanangle_x,accurateFit_[ievt].tanangle_y);
}

bool SignalRegion::fillEnergies(const unsigned ievt,
				const std::vector<HGCSSSamplingSection> & ssvec,
				const std::vector<HGCSSSimHit> & simhitvec,
				const std::vector<HGCSSRecoHit> & rechitvec,
				const unsigned nPuVtx){
  const FitResult & fit = accurateFit_[ievt];
  return fillEnergies(ievt,ssvec,simhitvec,rechitvec,nPuVtx,fit);
}

bool SignalRegion::fillEnergies(const unsigned ievt,
				const std::vector<HGCSSSamplingSection> & ssvec,
				const std::vector<HGCSSSimHit> & simhitvec,
				const std::vector<HGCSSRecoHit> & rechitvec,
				const unsigned nPuVtx,
				const FitResult & fit){
  
 
  //fill weights for first event only: same in all events
  if (firstEvent_){
    absweight_.clear();
    absweight_.reserve(nLayers_);
    std::cout << " -- Absorber weights used for total energy:" << std::endl;
    for(unsigned iL(0); iL<nLayers_; iL++){
      double w = ssvec[iL].volX0trans()/ssvec[1].volX0trans();
      std::cout << " - Layer " << iL << " w=" << w << std::endl;
      absweight_.push_back(w);
    }
    firstEvent_=false;
  }

  if (absweight_.size()!=nLayers_) {
    std::cout << " -- Error! Not all layers found! Only: " << absweight_.size() << ". Fix code." << std::endl;
    exit(1);
  }
  
    //initialise values for current event
  evtIdx_ = ievt;
  totalE_ = 0;
  wgttotalE_ = 0;
  
  for (unsigned iL(0); iL<nLayers_;++iL){
    for (unsigned iSR(0);iSR<nSR_;++iSR){
      energySR_[iL][iSR] = 0;
      subtractedenergySR_[iL][iSR] = 0;
    }
  }

  if(!fit.found) {
    //std::cout << " -- Event " << ievt << " skipped, accurate position not found." << std::endl;
    nSkipped_++;
    //fill tree to find correspondance between noPu and PU...
    outtree_->Fill();
    return false;
  }
  
  //std::cout << " -- Accurate direction for evt " << ievt << ": " << std::endl;
  //getAccurateDirection(ievt).Print();

  //initialise accuratepos per layer
  std::vector<ROOT::Math::XYZPoint> eventPos;
  eventPos.resize(nLayers_,ROOT::Math::XYZPoint(0,0,0));
  for (unsigned iL(0); iL<nLayers_;++iL){
    eventPos[iL] = getAccuratePos(fit,iL);
  }

  //double etacor = 1./tanh(getAccurateDirection(ievt).eta());

  //get event-by-event PU
  //get PU contrib from elsewhere in the event
  //loop over phi with same etamax
  //take average per layer: not all 9 cells of 3*3 area have hits...
  /*	std::vector<double> puE;
	puE.resize(nLayers_,0);
	if (nPuVtx>0){
	unsigned nRandomCones = 50;
	double phistep = TMath::Pi()/nRandomCones;
	for (unsigned ipm(0);ipm<nRandomCones;++ipm){
	std::vector<double> xmaxrc;
	xmaxrc.resize(nLayers_,0);
	std::vector<double> ymaxrc;
	ymaxrc.resize(nLayers_,0);
	double phirc = phimax-TMath::Pi();
	if (phirc < -1.*TMath::Pi()) phirc+=2.*TMath::Pi();
	if (ipm%2==0) phirc += ipm/2*phistep+phistep/2.;
	else  phirc = phirc - ipm/2*phistep-phistep/2.;
	if (phirc < -1.*TMath::Pi()) phirc+=2.*TMath::Pi();
	//take from geom to not be biased by hit having PU, because
	//not from geom means find cell with a hit closest to maxpos...
	getMaximumCellFromGeom(phirc,etamax,xmaxrc,ymaxrc);
	getPuContribution(rechitvec,xmaxrc,ymaxrc,puE);
	}
	
	//normalise to one cell: must count cells with 0 hit !
	//use cell size at etamax...
	for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
	unsigned nCells = nRandomCones*nSR_*geomConv_.cellSize()/geomConv_.cellSize(iL,etamax);
	puE[iL] = puE[iL]/nCells;
	}
	
	}//if PU
  */
  
  
  // Define different signal region and sum over energy
  for (unsigned iH(0); iH<rechitvec.size(); ++iH){//loop on hits
    const HGCSSRecoHit & lHit = rechitvec[iH];
    
    unsigned layer = lHit.layer();
    if (layer >= nLayers_) {
      continue;
    }
    double posx = lHit.get_x();
    if (fixForPuMixBug_) posx-=1.25;
    double posy = lHit.get_y();
    if (fixForPuMixBug_) posy-=1.25;
    double energy = lHit.energy();
    double leta = lHit.eta();

    //not interested in hits outside of acceptance...
    if (leta<1.4 || leta>3.0) continue;

    double etacor = fabs(1./tanh(leta));

    totalE_ += energy;
    wgttotalE_ += energy*absweight_[layer]*etacor;    
    
    double puE = puDensity_.getDensity(leta,layer,geomConv_.cellSizeInCm(layer,leta),nPuVtx);
    double subtractedenergy = std::max(0.,energy - puE);
    double halfCell = 0.5*geomConv_.cellSize(layer,leta);
    
    double dx = eventPos[layer].x()-posx;
    double dy = eventPos[layer].y()-posy;
    
    //SR0-4
    for (unsigned isr(0); isr<nSR_;++isr){
      if ( (fabs(dx) <= ((isr+1)*halfCell)) && (fabs(dy) <= ((isr+1)*halfCell))){
	energySR_[layer][isr] += energy*absweight_[layer]*etacor;
	subtractedenergySR_[layer][isr] += subtractedenergy*absweight_[layer]*etacor;
      }
    }
  }//loop on hits
  
  outputFile_->cd(outputDir_.c_str());
  fillHistograms();
  outtree_->Fill();
  return true;
}

void SignalRegion::finalise(){
  outputFile_->Flush();
  std::cout << " -- Histograms for signal regions have been filled !" << std::endl;
  std::cout << " -- Number of skipped events: " << nSkipped_ << std::endl;
}


void SignalRegion::initialiseHistograms(){

    outputFile_->cd(outputDir_.c_str());

    outtree_ = new TTree("Ereso","Tree to save energies in signal regions");


    outtree_->Branch("eventIndex",&evtIdx_);
    outtree_->Branch("rawEtotal",&totalE_);
    outtree_->Branch("wgtEtotal",&wgttotalE_);

    std::vector<double> emptyvec;
    emptyvec.resize(nSR_,0);
    energySR_.resize(nLayers_,emptyvec);
    subtractedenergySR_.resize(nLayers_,emptyvec);

    std::ostringstream label;
    for (unsigned iL(0); iL<nLayers_;++iL){
      for (unsigned iSR(0);iSR<nSR_;++iSR){
	label.str("");
	label << "energy_" << iL << "_SR" << iSR;
	outtree_->Branch(label.str().c_str(),&energySR_[iL][iSR]);
	label.str("");
	label << "subtractedenergy_" << iL << "_SR" << iSR;
	outtree_->Branch(label.str().c_str(),&subtractedenergySR_[iL][iSR]);
      }
    }

    //outtree_->SetCacheSize(100000000);
    //outtree_->SetCacheLearnEntries(100);
    //tree_ptr->LoadTree(evt);

    p_rawEtotal = new TH1F("p_rawEtotal", "Total E (MIP)", 5000,0,200000);
    p_wgtEtotal = new TH1F("p_wgtEtotal", "Total weighted E (MIP)",5000, 0, 200000);

    //p_rawESR.resize(nSR_,0);
    p_wgtESR.resize(nSR_,0);
    //p_rawSubtractESR.resize(nSR_,0);
    p_wgtSubtractESR.resize(nSR_,0);

    for (unsigned iSR(0);iSR<nSR_;++iSR){
      //label.str("");
      //label << "rawESR" << iSR;
      //p_rawESR[iSR] = new TH1F(("p_"+label.str()).c_str(),";E_{SR} (MIPs);events", 5000,0,200000);
      //p_rawESR[iSR]->StatOverflows();

      label.str("");
      label << "wgtESR" << iSR;
      p_wgtESR[iSR] = new TH1F(("p_"+label.str()).c_str(),";E_{SR} (MIPs);events", 5000,0,200000);
      p_wgtESR[iSR]->StatOverflows();
      
      //label.str("");
      //label << "rawSubtractESR" << iSR;
      //p_rawSubtractESR[iSR] = new TH1F(("p_"+label.str()).c_str(),";E_{SR}^{PUsubtr} (MIPs);events", 5000,0,200000);
      //p_rawSubtractESR[iSR]->StatOverflows();
      label.str("");
      label << "wgtSubtractESR" << iSR;
      p_wgtSubtractESR[iSR] = new TH1F(("p_"+label.str()).c_str(),";E_{SR}^{PUsubtr} (MIPs);events", 5000,0,200000);
      p_wgtSubtractESR[iSR]->StatOverflows();
    }//loop on sr

}

void SignalRegion::fillHistograms(){
        
  outputFile_->cd(outputDir_.c_str());

  p_rawEtotal->Fill(totalE_);
  p_wgtEtotal->Fill(wgttotalE_);

  for (unsigned iSR(0);iSR<nSR_;++iSR){
    //Fill energy without PU subtraction
    bool subtractPU = false;
    //p_rawESR[iSR]->Fill( getEtotalSR(iSR, subtractPU));

    double wgtESR = 0;
    for(unsigned iL(0); iL < nLayers_;iL++){
      wgtESR += getSR(iSR, iL, subtractPU);
    }
    p_wgtESR[iSR]->Fill(wgtESR);

    //Fill energy after PU subtraction
    subtractPU = true;
    //p_rawSubtractESR[iSR]->Fill( getEtotalSR(iSR, subtractPU));
    double wgtSubtractESR = 0;
    for(unsigned iL(0); iL < nLayers_;iL++){
      wgtSubtractESR += getSR(iSR, iL, subtractPU);
    }
    p_wgtSubtractESR[iSR]->Fill(wgtSubtractESR);

  }//loop on SR

}




