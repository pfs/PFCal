#include "EventAction.hh"

#include "RunAction.hh"
#include "EventActionMessenger.hh"
#include "DetectorConstruction.hh"

#include "HGCSSInfo.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

//
EventAction::EventAction()
{
  runAct = (RunAction*)G4RunManager::GetRunManager()->GetUserRunAction();
  eventMessenger = new EventActionMessenger(this);
  printModulo = 10;
  outF_=TFile::Open("PFcal.root","RECREATE");
  outF_->cd();

  double xysize = ((DetectorConstruction*)G4RunManager::GetRunManager()->GetUserDetectorConstruction())->GetCalorSizeXY();

  shape_ = ((DetectorConstruction*)G4RunManager::GetRunManager()->GetUserDetectorConstruction())->getShape();

  //save some info
  HGCSSInfo *info = new HGCSSInfo();
  info->calorSizeXY(xysize);
  info->cellSize(CELL_SIZE_X);
  info->model(((DetectorConstruction*)G4RunManager::GetRunManager()->GetUserDetectorConstruction())->getModel());
  info->version(((DetectorConstruction*)G4RunManager::GetRunManager()->GetUserDetectorConstruction())->getVersion());
  info->shape(shape_);
  std::cout << " -- check Info: version = " << info->version()
	    << " model = " << info->model()
	    << " shape = " << shape_
	    << std::endl;
  outF_->WriteObjectAny(info,"HGCSSInfo","Info");

  //honeycomb or diamond or triangles
  geomConv_ = new HGCSSGeometryConversion(info->model(),CELL_SIZE_X);
  if (shape_==2) geomConv_->initialiseDiamondMap(xysize,10.);
  else if (shape_==3) geomConv_->initialiseTriangleMap(xysize,10.*sqrt(2.));
  else if (shape_==1) geomConv_->initialiseHoneyComb(xysize,CELL_SIZE_X);
  else if (shape_==4) geomConv_->initialiseSquareMap(xysize,100.);
  //square map for FHCAL Scint + BH Scint
  double etamin = ((DetectorConstruction*)G4RunManager::GetRunManager()->GetUserDetectorConstruction())->GetMinEta();
  double etamax = ((DetectorConstruction*)G4RunManager::GetRunManager()->GetUserDetectorConstruction())->GetMaxEta();
  std::vector<G4double> phiSegmentationList=((DetectorConstruction*)G4RunManager::GetRunManager()->GetUserDetectorConstruction())->GetPhiSectorSegmentation();
  geomConv_->initialiseSquareMap1(etamin,etamax,-1.*TMath::Pi(),TMath::Pi(),phiSegmentationList[0]);//eta phi segmentation
  geomConv_->initialiseSquareMap2(etamin,etamax,-1.*TMath::Pi(),TMath::Pi(),phiSegmentationList[1]);
  
  //for the TDR geometry this is the first layer to use a coarser eta-phi segmentation
  SetFirstBHCAL2Layer(((DetectorConstruction*)G4RunManager::GetRunManager()->GetUserDetectorConstruction())->GetFirstBHCAL2Layer());
  std::cout << "[EventAction] first BHCAL2 layer will be set to " << GetFirstBHCAL2Layer() << std::endl
            << "              phi segmentation in BHCAL1 set to " << phiSegmentationList[0] << std::endl
            << "              phi segmentation in BHCAL2 set to " << phiSegmentationList[1] << std::endl;

  tree_=new TTree("HGCSSTree","HGC Standalone simulation tree");
  tree_->Branch("HGCSSEvent","HGCSSEvent",&event_);
  tree_->Branch("HGCSSSamplingSectionVec","std::vector<HGCSSSamplingSection>",&ssvec_);
  tree_->Branch("HGCSSSimHitVec","std::vector<HGCSSSimHit>",&hitvec_);
  tree_->Branch("HGCSSAluSimHitVec","std::vector<HGCSSSimHit>",&alhitvec_);
  tree_->Branch("HGCSSGenParticleVec","std::vector<HGCSSGenParticle>",&genvec_);

  //fout_.open("ProcessDepAbove5MeV.dat");
  //if (!fout_.is_open()){
  //std::cout << " -- Output file could not be opened..." << std::endl;
  //exit(1);
  //}
}

//
EventAction::~EventAction()
{
  outF_->cd();
  tree_->Write();
  outF_->Close();
  //fout_.close();
  delete eventMessenger;
}

//
void EventAction::BeginOfEventAction(const G4Event* evt)
{  
  evtNb_ = evt->GetEventID();
  if (evtNb_%printModulo == 0) { 
    G4cout << "\n---> Begin of event: " << evtNb_ << G4endl;
    CLHEP::HepRandom::showEngineStatus();
  }
  //fout_ << "Event " << evtNb_ << std::endl;

}

//
void EventAction::Detect(G4double edep, G4double stepl,G4double globalTime, 
			 G4int pdgId, G4VPhysicalVolume *volume, const G4ThreeVector & position, 
			 G4int trackID, G4int parentID,
			 const HGCSSGenParticle & genPart)
{
  for(size_t i=0; i<detector_->size(); i++) (*detector_)[i].add(edep,stepl,globalTime,pdgId,volume,position,trackID,parentID,i);
  if (genPart.isIncoming()) genvec_.push_back(genPart);
}

bool EventAction::isFirstVolume(const std::string volname) const{
  if (detector_->size()>0 && (*detector_)[0].n_elements>0){
    bool found = false;
    //for (unsigned iS(0); iS<(*detector_)[0].n_sectors;++iS){
      //if ((((*detector_)[0].ele_vol[(*detector_)[0].n_elements*iS])->GetName())==volname.c_str() || (*detector_)[0].supportcone_vol->GetName()==volname.c_str()) found = true;
      //}
    if ((*detector_)[0].dummylayer_vol->GetName()==volname.c_str()) found = true;
    return found;
  }
  return "";
}

//
void EventAction::EndOfEventAction(const G4Event* g4evt)
{
  //return;
  bool debug(evtNb_%printModulo == 0);
  hitvec_.clear();

  event_.eventNumber(evtNb_);

  //std::cout << " -- Number of primary vertices: " << g4evt->GetNumberOfPrimaryVertex() << std::endl
  //<< " -- vtx pos x=" << g4evt->GetPrimaryVertex(0)->GetX0() 
  //	    << " y=" << g4evt->GetPrimaryVertex(0)->GetY0()
  //	    << " z=" << g4evt->GetPrimaryVertex(0)->GetZ0()
  //	    << " t=" << g4evt->GetPrimaryVertex(0)->GetT0()
  //	    << std::endl;

  event_.vtx_x(g4evt->GetPrimaryVertex(0)->GetX0());
  event_.vtx_y(g4evt->GetPrimaryVertex(0)->GetY0());
  event_.vtx_z(g4evt->GetPrimaryVertex(0)->GetZ0());
  //event_.cellSize(CELL_SIZE_X);

  ssvec_.clear();
  ssvec_.reserve(detector_->size());

  for(size_t i=0; i<detector_->size(); i++) 
    {
      HGCSSSamplingSection lSec;
      lSec.volNb(i);
      lSec.sensitiveZ((*detector_)[i].sensitiveZ);
      lSec.volX0trans((*detector_)[i].getAbsorberX0());
      lSec.voldEdx((*detector_)[i].getAbsorberdEdx());
      lSec.volLambdatrans((*detector_)[i].getAbsorberLambda());
      lSec.absorberE((*detector_)[i].getAbsorbedEnergy());
      lSec.measuredE((*detector_)[i].getMeasuredEnergy(false));
      lSec.totalE((*detector_)[i].getTotalEnergy());
      lSec.gFrac((*detector_)[i].getPhotonFraction());
      lSec.eFrac((*detector_)[i].getElectronFraction());
      lSec.muFrac((*detector_)[i].getMuonFraction());
      lSec.neutronFrac((*detector_)[i].getNeutronFraction());
      lSec.hadFrac((*detector_)[i].getHadronicFraction());
      lSec.avgTime((*detector_)[i].getAverageTime());
      lSec.nSiHits((*detector_)[i].getTotalSensHits());
      ssvec_.push_back(lSec);
      if (evtNb_==1) std::cout << "if (layer==" << i << ") return " 
			       <<  lSec.voldEdx() << ";"
			       << std::endl;
      //std::cout << " n_sens_ele = " << (*detector_)[i].n_sens_elements << std::endl;
      bool is_scint = (*detector_)[i].hasScintillator;
      for (unsigned idx(0); idx<(*detector_)[i].n_sens_elements; ++idx){
	std::map<unsigned,HGCSSSimHit> lHitMap;
	std::pair<std::map<unsigned,HGCSSSimHit>::iterator,bool> isInserted;
	
	//if (i>0) (*detector_)[i].trackParticleHistory(idx,(*detector_)[i-1].getSiHitVec(idx));
	
	//std::cout << " si layer " << idx << " " << (*detector_)[i].getSiHitVec(idx).size() << std::endl;

	for (unsigned iSiHit(0); iSiHit<(*detector_)[i].getSiHitVec(idx).size();++iSiHit){
	  G4SiHit lSiHit = (*detector_)[i].getSiHitVec(idx)[iSiHit];
	  HGCSSSimHit lHit(lSiHit,idx,is_scint? (i<firstBHCAL2Layer_?geomConv_->squareMap1():geomConv_->squareMap2()): (shape_==4 ?geomConv_->squareMap() : shape_==2?geomConv_->diamondMap():shape_==3?geomConv_->triangleMap():geomConv_->hexagonMap()),CELL_SIZE_X,is_scint?true:false);
	  
	  isInserted = lHitMap.insert(std::pair<unsigned,HGCSSSimHit>(lHit.cellid(),lHit));
	  if (!isInserted.second) isInserted.first->second.Add(lSiHit);
	}
	std::map<unsigned,HGCSSSimHit>::iterator lIter = lHitMap.begin();
	hitvec_.reserve(hitvec_.size()+lHitMap.size());
	for (; lIter != lHitMap.end(); ++lIter){
	  (lIter->second).calculateTime();
	  hitvec_.push_back(lIter->second);
	}

      }//loop on sensitive layers

      //support cone
      std::map<unsigned,HGCSSSimHit> lHitMap;
      std::pair<std::map<unsigned,HGCSSSimHit>::iterator,bool> isInserted;
      for (unsigned iAlHit(0); iAlHit<(*detector_)[i].getAlHitVec().size();++iAlHit){
	G4SiHit lAlHit = (*detector_)[i].getAlHitVec()[iAlHit];
	HGCSSSimHit lHit(lAlHit,0,geomConv_->squareMap());
	isInserted = lHitMap.insert(std::pair<unsigned,HGCSSSimHit>(0,lHit));
	if (!isInserted.second) isInserted.first->second.Add(lAlHit);
      }
      std::map<unsigned,HGCSSSimHit>::iterator lIter = lHitMap.begin();
      alhitvec_.reserve(alhitvec_.size()+lHitMap.size());
      for (; lIter != lHitMap.end(); ++lIter){
	(lIter->second).calculateTime();
	alhitvec_.push_back(lIter->second);
      }

      if(debug) {
	(*detector_)[i].report( (i==0) );
      }
      //if (i==0) G4cout << " ** evt " << evt->GetEventID() << G4endl;
      (*detector_)[i].resetCounters();

    }
  if(debug){
    G4cout << " -- Number of truth particles = " << genvec_.size() << G4endl
	   << " -- Number of simhits = " << hitvec_.size() << G4endl
	   << " -- Number of aluminium simhits = " << alhitvec_.size() << G4endl
	   << " -- Number of sampling sections = " << ssvec_.size() << G4endl;
    
  }

  tree_->Fill();
  
  //reset vectors
  genvec_.clear();
  hitvec_.clear();
  alhitvec_.clear();
  ssvec_.clear();
}
