#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<map>
#include <boost/algorithm/string.hpp>
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "boost/function.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"

#include "HGCSSEvent.hh"
#include "HGCSSInfo.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSParameters.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDigitisation.hh"
#include "HGCSSDetector.hh"
#include "HGCSSGeometryConversion.hh"
#include "HGCSSPUenergy.hh"
#include "HGCSSMipHit.hh"

#include "PositionFit.hh"
#include "SignalRegion.hh"

#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"

#include "utilities.h"

using boost::lexical_cast;
namespace po=boost::program_options;

int main(int argc, char** argv){//main  

  //Input output and config options
  std::string cfg;
  unsigned pNevts;
  std::string inFilePath;
  std::string outFilePath;
  unsigned nRuns;
  unsigned startIdx;
  unsigned debug;
  double etamean;
  double deta;
  po::options_description preconfig("Configuration"); 
  preconfig.add_options()("cfg,c",po::value<std::string>(&cfg)->required());
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(preconfig).allow_unregistered().run(), vm);
  po::notify(vm);
  po::options_description config("Configuration");
  config.add_options()
    //Input output and config options //->required()
    ("pNevts,n",       po::value<unsigned>(&pNevts)->default_value(0))
    ("inFilePath,i",   po::value<std::string>(&inFilePath)->required())
    ("outFilePath,o",  po::value<std::string>(&outFilePath)->required())
    ("nRuns",        po::value<unsigned>(&nRuns)->default_value(0))
    ("startIdx",     po::value<unsigned>(&startIdx)->default_value(0))
    ("debug,d",        po::value<unsigned>(&debug)->default_value(0))
    ("etamean,e",      po::value<double>(&etamean)->default_value(2.8))
    ("deta",      po::value<double>(&deta)->default_value(0.05))
    ;

  po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
  po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
  po::notify(vm);


  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << inFilePath << std::endl
	    << " -- Output file path: " << outFilePath << std::endl
	    << " -- mean eta: " << etamean 
	    << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events per run." << std::endl;

  TRandom3 lRndm(0);
  std::cout << " -- Random number seed: " << lRndm.GetSeed() << std::endl;


  /////////////////////////////////////////////////////////////
  //hardcoded
  /////////////////////////////////////////////////////////////

  //fix for eta=2.8 digi which did not have noise only added...
  const bool addNoiseOnly = false;//etamean>=2.7?true:false;

  //to run faster....
  const bool filter0noise = true;
  //works for both small and large cells...

  const double neighRadius = etamean >= 2.4 ? 9 : 13;//mm for 6 closest neighbours

  const double sizeFH = 2./360.;//0.01745;
  const double sizeBH = 2./288.;//0.02182;

  //global threshold to reduce size of noise hits
  const double threshMin = 0.5;
  //save only central hits with 0.5 mips...
  const double threshCentralHitMin = 0.5;
  const double threshCentralHitMax = 5;
  const double threshBeforeAfterMin = 0.5;
  const double threshBeforeAfterMax = 5;
  const double threshCentralNeighMax = 3;//3*sigma_noise...

  std::cout << " ---- Selection settings: ---- " << std::endl
	    << " -------threshMin " << threshMin << std::endl
	    << " -------threshCentralHitMin " << threshCentralHitMin << std::endl
	    << " -------threshCentralHitMax " << threshCentralHitMax << std::endl
	    << " -------threshBeforeAfterMin " << threshBeforeAfterMin << std::endl
	    << " -------threshBeforeAfterMax " << threshBeforeAfterMax << std::endl
	    << " -------threshCentralNeighMax " << threshCentralNeighMax << " * sigma_Noise" << std::endl
	    << " ------------------------------------------" << std::endl;


  //const unsigned nEta = 1;

  const unsigned nNoise = 8;//5;//10;
  const double noise[nNoise] = {0,0.13,0.2,0.27,0.35,0.4,0.45,0.5};
  //double eta[nEta];// = {2.85};//1.7,2.0,2.5};

  std::cout << "Noise values: " ;
  for (unsigned iN(0); iN<nNoise;++iN){
    std::cout << noise[iN] << " ";
  }
  std::cout << std::endl;

  //eta[0] = etamean;

  //const double deta = etamean>=2.7?0.05:0.02;

  //const double Ethresh = 0.6;

  /////////////////////////////////////////////////////////////
  //input
  /////////////////////////////////////////////////////////////
  TChain *lTree = new TChain("RecoTree");
  TFile * recFile = 0;
  TFile * lastrecFile = 0;
  if (nRuns == 0){
    if (!testInputFile(inFilePath,recFile)) return 1;
    lTree->AddFile(inFilePath.c_str());
    lastrecFile = recFile;
  }
  else {
    for (unsigned i(startIdx);i<startIdx+nRuns;++i){
      std::ostringstream lstrrec;
      lstrrec << inFilePath << "_run" << i << ".root";
      if (!testInputFile(lstrrec.str(),recFile)) continue;
      lTree->AddFile(lstrrec.str().c_str());
      lastrecFile = recFile;
    }
  }
  if (!lTree){
    std::cout << " -- Error, tree RecoTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  std::cout << " Trees added." << std::endl;


  HGCSSInfo * info=(HGCSSInfo*)lastrecFile->Get("Info");

  assert(info);
  //info->Print();

  const unsigned versionNumber = info->version();
  const unsigned model = info->model();
  //CAMM-dirty fix
  //const unsigned shape = inputStr.find("Diamond")!=inputStr.npos? 2: inputStr.find("Triangle")!=inputStr.npos?3 :1 ;//info->shape();
  const unsigned shape = info->shape();
  //const double cellSize = 6.496345;//info->cellSize();
  //const double calorSizeXY = 2800*2;//info->calorSizeXY();
  const double cellSize = info->cellSize();
  const double calorSizeXY = info->calorSizeXY();

  bool isTBsetup = (model != 2);
  bool bypassR = false;
  if (isTBsetup) bypassR = true;

  HGCSSDetector & myDetector = theDetector();
  myDetector.buildDetector(versionNumber,true,false,bypassR);

  //corrected for Si-Scint overlap
  const unsigned nLayers = 52;//etamean<2.3? myDetector.nLayers(): 52;


  std::cout << " -- Calor size XY = " << calorSizeXY
	    << ", version number = " << versionNumber 
	    << ", model = " << model << std::endl
	    << " -- cellSize = " << cellSize
	    << ", shape = " << shape
	    << ", nLayers = " << nLayers
	    << std::endl;
  HGCSSGeometryConversion geomConv(model,cellSize,bypassR,3);

  geomConv.setXYwidth(calorSizeXY);
  geomConv.setVersion(versionNumber);
  //const double xWidth = geomConv.getXYwidth();
  
  if (shape==2) geomConv.initialiseDiamondMap(calorSizeXY,10.);
  else if (shape==3) geomConv.initialiseTriangleMap(calorSizeXY,10.*sqrt(2.));
  else if (shape==1) geomConv.initialiseHoneyComb(calorSizeXY,cellSize);
  else if (shape==4) geomConv.initialiseSquareMap(calorSizeXY,10.);

  //square map for BHCAL
  //square map for BHCAL
  geomConv.initialiseSquareMap1(1.3,3.0,-1.*TMath::Pi(),1.*TMath::Pi(),sizeFH);//eta phi segmentation
  geomConv.initialiseSquareMap2(1.3,3.0,-1.*TMath::Pi(),1.*TMath::Pi(),sizeBH);//eta phi segmentation
  std::vector<unsigned> granularity;
  granularity.resize(myDetector.nLayers(),1);
  geomConv.setGranularity(granularity);
  geomConv.initialiseHistos();

  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  ///////// Output File // /////////////////////////
  //////////////////////////////////////////////////
  //////////////////////////////////////////////////

  TFile *outputFile = TFile::Open(outFilePath.c_str(),"RECREATE");
  
  if (!outputFile) {
    std::cout << " -- Error, output file " << outFilePath << " cannot be opened. Please create output directory. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
  }
  outputFile->cd();

  TTree *miptree=new TTree("MipTree","HGC MipAnalysis tree");

  HGCSSGenParticleVec lGenParts;
  miptree->Branch("HGCSSGenParticleVec","std::vector<HGCSSGenParticle>",&lGenParts);


  std::vector<HGCSSMipHitVec> miphitvec;
  HGCSSMipHitVec init;
  //need to init such that the vector is not moved around in memory!!
  HGCSSMipHit dummy;
  const unsigned nHitsInit = 20000;
  init.reserve(nHitsInit);//,dummy);
  miphitvec.resize(nNoise,init);

  for (unsigned iN(0); iN<nNoise;++iN){
    std::ostringstream label;
    //root doesn't like . in branch names.....
    label << "MipHitVec_noise" << iN;//noise[iN];
    //std::cout << " vec pointer branch " << iN << " " << &miphitvec[iN] << std::endl;
     miptree->Branch(label.str().c_str(),"std::vector<HGCSSMipHit>",&miphitvec[iN]);
  }

  const unsigned nEvts = ((pNevts > lTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lTree->GetEntries()) : pNevts) ;
  
  std::cout << " -- Processing " << nEvts << " events out of " << lTree->GetEntries() << std::endl;


  //loop on events
  HGCSSEvent * event = 0;
  std::vector<HGCSSRecoHit> * rechitvec = 0;
  std::vector<HGCSSGenParticle> * genvec = 0;
  unsigned nPuVtx = 0;

  lTree->SetBranchAddress("HGCSSEvent",&event);
  lTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  if (lTree->GetBranch("HGCSSGenParticleVec")) lTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);
  else genvec = new std::vector<HGCSSGenParticle>;
  if (lTree->GetBranch("nPuVtx")) lTree->SetBranchAddress("nPuVtx",&nPuVtx);

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
    lTree->GetEntry(ievt);

    //save genparticles
    lGenParts.clear();
    lGenParts.reserve((*genvec).size());
    for (unsigned iP(0); iP<(*genvec).size(); ++iP){//loop on gen particles
      HGCSSGenParticle lPart = (*genvec)[iP];
      lGenParts.push_back(lPart);
    }


    //add noise only contributions: need to be done consistently for the event... i.e. same cell can be neighbour to several hits....
    if (addNoiseOnly){
      //fill map
      bool isScint = false;
      if (debug) std::cout << " - Event contains " << (*rechitvec).size() << " rechits." << std::endl;
      unsigned notInEtaRange = 0;
      unsigned notInLayerRange = 0;
      for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
	HGCSSRecoHit lHit = (*rechitvec)[iH];
	double leta = lHit.eta();
	if (debug>1) std::cout << " -- hit " << iH << " eta " << leta << std::endl; 
	//clean up rechit collection
	if (fabs(leta-etamean)>= deta){
	  rechitvec->erase(rechitvec->begin()+iH);
	  --iH;
	  notInEtaRange++;
	  continue;
	}
	unsigned layer = lHit.layer();
	if (layer>=myDetector.nLayers()) {
	  rechitvec->erase(rechitvec->begin()+iH);
	  --iH;
	  notInLayerRange++;
	  continue;
	}
	const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
	isScint = subdet.isScint;
	TH2Poly *map = isScint?(subdet.type==DetectorEnum::BHCAL1?geomConv.squareMap1():geomConv.squareMap2()): shape==4?geomConv.squareMap() : shape==2?geomConv.diamondMap() : shape==3? geomConv.triangleMap(): geomConv.hexagonMap();

	unsigned cellid = 0;
	if (isScint){
	  ROOT::Math::XYZPoint pos = ROOT::Math::XYZPoint(lHit.get_x(),lHit.get_y(),lHit.get_z());
	  cellid = map->FindBin(pos.eta(),pos.phi());
	} else {
	  cellid = map->FindBin(lHit.get_x(),lHit.get_y());
	}
	geomConv.fill(lHit.layer(),lHit.energy(),0,cellid,lHit.get_z());
      }//loop on hits

      if (debug) std::cout << " - In eta range, event contains " << (*rechitvec).size() << " rechits." << std::endl;

      unsigned nTotBins = 0;
      unsigned nTotCheck = 0;
      unsigned nAddedCheck = 0;
      for (unsigned iL(0); iL<myDetector.nLayers(); ++iL){//loop on layers
	std::map<unsigned,MergeCells> & histE = geomConv.get2DHist(iL);
	//std::cout << " Check layer " << iL << " nHits = " << histE.size() << std::endl;
	nTotCheck += histE.size();
	const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(iL);
	bool isScint = subdet.isScint;

	std::map<int,std::pair<double,double> > & geom = isScint?(subdet.type==DetectorEnum::BHCAL1?geomConv.squareGeom1:geomConv.squareGeom2): shape==4?geomConv.squareGeom:shape==2?geomConv.diamGeom:shape==3?geomConv.triangleGeom:geomConv.hexaGeom;

	unsigned nBins = geom.size();
	nTotBins += nBins;
	double meanZpos = myDetector.sensitiveZ(iL);

	//extend map to include all cells in eta ring
	for (unsigned iB(1); iB<nBins+1;++iB){
	  std::pair<double,double> xy = geom[iB];
	  if (isScint) {
	    HGCSSGeometryConversion::convertFromEtaPhi(xy,meanZpos);
	  }
	  ROOT::Math::XYZPoint lpos = ROOT::Math::XYZPoint(xy.first,xy.second,meanZpos);
	  double eta = lpos.eta();
	  bool passeta = fabs(eta-etamean)<deta;
	  if (!passeta) continue;
	  MergeCells tmpCell;
	  tmpCell.energy = 0;
	  tmpCell.time = 0;
	  std::pair<std::map<unsigned,MergeCells>::iterator,bool> isInserted = histE.insert(std::pair<unsigned,MergeCells>(iB,tmpCell));
	  if (isInserted.second==true){
	    HGCSSRecoHit lRecHit;
	    lRecHit.layer(iL);
	    lRecHit.energy(0);
	    lRecHit.adcCounts(0);
	    lRecHit.x(xy.first);
	    lRecHit.y(xy.second);
	    lRecHit.z(meanZpos);
	    lRecHit.noiseFraction(1);
	    rechitvec->push_back(lRecHit);
	    nAddedCheck++;
	  }
	}
      }//loop on layers
      if (debug) {
	std::cout << " - Adding noise hits, event now contains " << (*rechitvec).size() << " rechits, not in eta range: " << notInEtaRange 
		  << " not in layer range (0-" << myDetector.nLayers() << "): " << notInLayerRange
		  << std::endl;
	//std::cout << " nSignalCheck = " << nTotCheck << " nAddedCheck = " << nAddedCheck << " " << nTotCheck+nAddedCheck << std::endl;
      }
    }//addNoiseOnly


    //loop once more to add noise and apply threshold to speed up processing
    //will remove half noise-only hits with E<0...

    std::vector<double> hitEnergy[nNoise];
    for (unsigned iN(0); iN<nNoise;++iN){
      miphitvec[iN].clear();
      miphitvec[iN].reserve(nHitsInit);
      //std::cout << " vec pointer event loop " << iN << " " << &miphitvec[iN] << " size " <<  miphitvec[iN].size() << std::endl;
      hitEnergy[iN].clear();
      hitEnergy[iN].resize((*rechitvec).size(),-1);
    }
    
    if (debug) std::cout << " - Event contains " << (*rechitvec).size() << " rechits." << std::endl;

    double zpos[myDetector.nLayers()];
    unsigned newLid[myDetector.nLayers()];
    for (unsigned iL(0); iL<myDetector.nLayers(); ++iL){//loop on layers
      zpos[iL] = myDetector.sensitiveZ(iL);
      if (iL<53) newLid[iL] = iL;
      else newLid[iL] = iL-17;
    }

    std::vector<bool> isSelected;
    isSelected.resize((*rechitvec).size(),false);
    unsigned nsel = 0;
    unsigned notInEtaRange = 0;
    unsigned notInLayerRange = 0;
    for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
      HGCSSRecoHit lHit = (*rechitvec)[iH];
      double leta = lHit.eta();
      //@FIXME should be etaphys....
      if (fabs(leta-etamean)>= deta){
	rechitvec->erase(rechitvec->begin()+iH);
	--iH;
	notInEtaRange++;
	//cut already applied at digi, should be in etaphys....
	continue;
      }
      //unsigned layer = lHit.layer();
      //if (layer>=myDetector.nLayers()) {
	//rechitvec->erase(rechitvec->begin()+iH);
	//--iH;
	//notInLayerRange++;
	//continue;
      //}

      
      if (filter0noise && lHit.energy()<0.1) {
	rechitvec->erase(rechitvec->begin()+iH);
	--iH;
	continue;
      }
      bool passOne = false;
      for (unsigned iN(0); iN<nNoise;++iN){
	double lNoise = lRndm.Gaus(0,noise[iN]);
	
	hitEnergy[iN][iH] = lHit.energy()+lNoise;
	if (hitEnergy[iN][iH]>threshMin) passOne=true;
      }
      if (!passOne) {
	rechitvec->erase(rechitvec->begin()+iH);
	--iH;
	continue;
      }
      isSelected[iH] = true;
      nsel++;
      if (debug && nsel%1000==0) std::cout << " -- hit " << iH/1000 << "k / " << (*rechitvec).size() << std::endl;
    }//loop on hits

    if (debug) {
      if (!addNoiseOnly) std::cout << " - not in eta range: " << notInEtaRange << std::endl
				   << " - not in layer range (0-" << myDetector.nLayers() << "): " << notInLayerRange
				   << std::endl;
      std::cout << " - After eta and minE threshold event contains " << (*rechitvec).size() << " rechits, " << nsel << " selected." << std::endl;
    }
    //if ((*rechitvec).size()>20000) {
    //geomConv.initialiseHistos();
    //continue;
    //}

    unsigned passSel[3] = {0,0,0};
    for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
      if (!filter0noise && !isSelected[iH]) continue;
      if (debug>1) std::cout << "...... Processing hit: " << iH << std::endl;
      else if (debug==1 && iH%1000 == 0) std::cout << "...... Processing hit: " << iH/1000 << "k" << std::endl;
      bool pass[nNoise];
      HGCSSRecoHit lCentralHit = (*rechitvec)[iH];
      bool passOne = false;
      for (unsigned iN(0); iN<nNoise;++iN){
	pass[iN] = true;
	if (hitEnergy[iN][iH]<threshCentralHitMin || hitEnergy[iN][iH]>threshCentralHitMax) {
	  pass[iN]=false;	  
	}
	else passOne=true;
	if (iN==0 && pass[0]) passSel[0]++;
      }
      if (!passOne) continue;

      unsigned hitlayer = lCentralHit.layer();
      unsigned newhitlayer = newLid[lCentralHit.layer()];
      double hitposx = lCentralHit.get_x();
      double hitposy = lCentralHit.get_y();
      double hitposz = lCentralHit.get_z();
      double etaPhys = getEtaPhys(hitposx,hitposy,hitposz,event);

      double hiteta = etaPhys;//lCentralHit.eta();
      double hitphi = lCentralHit.phi();
      double hittheta = 2*atan(exp(-hiteta));
      double sincos = sin(hittheta)*cos(hitphi);
      double sinsin = sin(hittheta)*sin(hitphi);
      HGCSSMipHit myHit[nNoise];
      for (unsigned iN(0); iN<nNoise;++iN){
	myHit[iN].setE(hitEnergy[iN][iH]);
	myHit[iN].setx(hitposx);
	myHit[iN].sety(hitposy);
	myHit[iN].setz(hitposz);
	myHit[iN].setLayer(hitlayer);
	myHit[iN].setnoiseFrac(lCentralHit.noiseFraction());
      }
      //find ref position of hit with same eta in previous and next layers
      std::pair<double,double> xyRef[5];
      for (int il(-2); il<3;++il){
	int neighlay = newhitlayer+il;
	if (neighlay<0) neighlay+=5;
	if (neighlay>(nLayers-1)) neighlay-=5;
	//bridge for the si-scint divide
	int oldneighlay = hitlayer<52? neighlay : hitlayer+il;
	if (oldneighlay>(myDetector.nLayers()-1)) oldneighlay-=5;
	if (hitlayer>52 && oldneighlay<53) oldneighlay -= 17;
	double rRef = zpos[oldneighlay]/cos(hittheta);
	double xRef = rRef*sincos;
	double yRef = rRef*sinsin;
	const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(oldneighlay);
	bool isScint = subdet.isScint;
	TH2Poly *map = isScint?(subdet.type==DetectorEnum::BHCAL1?geomConv.squareMap1():geomConv.squareMap2()): shape==4?geomConv.squareMap(): geomConv.hexagonMap();
	std::map<int,std::pair<double,double> > & geom = isScint?(subdet.type==DetectorEnum::BHCAL1?geomConv.squareGeom1:geomConv.squareGeom2):geomConv.hexaGeom;
	unsigned cellid = isScint?map->FindBin(hiteta,hitphi):map->FindBin(xRef,yRef);
	xyRef[il+2] = geom[cellid];
	myHit[0].set_xref_neighlay(il+2,xyRef[il+2].first);
	myHit[0].set_yref_neighlay(il+2,xyRef[il+2].second);
      }

      //loop on hits to find neighbours
      for (unsigned jH(0); jH<(*rechitvec).size(); ++jH){//loop on hits
	if (!filter0noise && !isSelected[jH]) continue;
	if (jH==iH) continue;
	HGCSSRecoHit lHit = (*rechitvec)[jH];
	unsigned newlayer = newLid[lHit.layer()];
	if ( (newhitlayer==0 && abs(newlayer-newhitlayer)<5)||
	     (newhitlayer==1 && abs(newlayer-newhitlayer)<4) ||
	     (newhitlayer==nLayers-1 && abs(newlayer-newhitlayer)<5) || 
	     (newhitlayer==nLayers-2 && abs(newlayer-newhitlayer)<4) || 
	     (newhitlayer>1 && newhitlayer<(nLayers-2) && abs(newlayer-newhitlayer)<3)
	     ){

	  const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(lHit.layer());
	  bool isScint = subdet.isScint;
	  //multiply by 1.1 to avoid rounding issues
	  double sizeEtaPhi = subdet.type==DetectorEnum::BHCAL1?sizeFH*1.1 : sizeBH*1.1;
	  //double rxyz = sqrt(pow(lHit.get_x(),2)+pow(lHit.get_y(),2)+pow(lHit.get_z(),2));
	  int lidx = newlayer-newhitlayer+2;
	  if (lidx>4) lidx-=5;
	  else if (lidx<0) lidx+=5;
	  double dx = isScint? lHit.eta()-xyRef[lidx].first : lHit.get_x()-xyRef[lidx].first;
	  double dy = isScint? lHit.phi()-xyRef[lidx].second : lHit.get_y()-xyRef[lidx].second;
	  //double dx = lHit.get_x()-hitposx;
	  //double dy = lHit.get_y()-hitposy;
	  double radius = sqrt(pow(dx,2)+pow(dy,2));
	  if ((isScint && fabs(dx)<sizeEtaPhi && fabs(dy)<sizeEtaPhi) || (!isScint && radius<neighRadius)){
	    //use -0.1 to protect against double precision
	    int index = 0;
	    if (isScint){
	      if (dy<-0.01) index = dx>0.01? 3 : dx<-0.01 ? 1 : 2;
	      else if (dy>0.01) index = dx>0.01? 8 : dx<-0.01 ? 6 : 7;
	      else index = dx>0.01? 5 : dx<-0.01 ? 4 : 0;
	    } else {
	      index = dy<-0.1? (dx<-0.1?1:fabs(dx)<0.01?2:3) : (dx<-0.1?4:fabs(dx)<0.01?5:6);
	    }

	    if (fabs(dx)<0.01&&fabs(dy)<0.01) index = 0;
	    if (debug && isScint && iH%1000==0){
	      std::cout << " Hit " << iH << " l " << hitlayer 
			<< " nl " << lHit.layer()
			<< " dx " << dx << " dy " << dy 
			<< " r " << radius << " idx " << index
			<< " xy " << lHit.get_x() << "," << lHit.get_y()
			<< " xyref " << xyRef[lidx].first << "," << xyRef[lidx].second
			<< " eta,phi " << lHit.eta() << "," << lHit.phi()
			<< " eta,phiref " << hiteta << "," << hitphi
			<< std::endl;
	    }
	    for (unsigned iN(0); iN<nNoise;++iN){
	      if (newlayer==newhitlayer) {
		//fill neighbour energies
		myHit[iN].set_neigh_e_samelayer(index,hitEnergy[iN][jH]);
	      }
	      if (newhitlayer==0){
		if (newlayer == (newhitlayer+2)) myHit[iN].set_neigh_e_prevlayer(index,hitEnergy[iN][jH]);
		if (newlayer == (newhitlayer+1)) myHit[iN].set_neigh_e_nextlayer(index,hitEnergy[iN][jH]);
		if (newlayer == (newhitlayer+3)) myHit[iN].set_neigh_e_prev2layer(index,hitEnergy[iN][jH]);
		if (newlayer == (newhitlayer+4)) myHit[iN].set_neigh_e_next2layer(index,hitEnergy[iN][jH]);
	      } else if (newhitlayer==1) {
		if (newlayer == (newhitlayer-1)) myHit[iN].set_neigh_e_prevlayer(index,hitEnergy[iN][jH]);
		if (newlayer == (newhitlayer+1)) myHit[iN].set_neigh_e_nextlayer(index,hitEnergy[iN][jH]);
		if (newlayer == (newhitlayer+2)) myHit[iN].set_neigh_e_prev2layer(index,hitEnergy[iN][jH]);
		if (newlayer == (newhitlayer+3)) myHit[iN].set_neigh_e_next2layer(index,hitEnergy[iN][jH]);
	      } else if (newhitlayer==nLayers-1){
		if (newlayer == (newhitlayer-1)) myHit[iN].set_neigh_e_prevlayer(index,hitEnergy[iN][jH]);
		if (newlayer == (newhitlayer-2)) myHit[iN].set_neigh_e_nextlayer(index,hitEnergy[iN][jH]);
		if (newlayer == (newhitlayer-3)) myHit[iN].set_neigh_e_prev2layer(index,hitEnergy[iN][jH]);
		if (newlayer == (newhitlayer-4)) myHit[iN].set_neigh_e_next2layer(index,hitEnergy[iN][jH]);
	      } else if (newhitlayer==nLayers-2) {
		if (newlayer == (newhitlayer-1)) myHit[iN].set_neigh_e_prevlayer(index,hitEnergy[iN][jH]);
		if (newlayer == (newhitlayer+1)) myHit[iN].set_neigh_e_nextlayer(index,hitEnergy[iN][jH]);
		if (newlayer == (newhitlayer-2)) myHit[iN].set_neigh_e_prev2layer(index,hitEnergy[iN][jH]);
		if (newlayer == (newhitlayer-3)) myHit[iN].set_neigh_e_next2layer(index,hitEnergy[iN][jH]);
	      } else {
		if (newlayer == (newhitlayer-1)) myHit[iN].set_neigh_e_prevlayer(index,hitEnergy[iN][jH]);
		if (newlayer == (newhitlayer+1)) myHit[iN].set_neigh_e_nextlayer(index,hitEnergy[iN][jH]);
		if (newlayer == (newhitlayer-2)) myHit[iN].set_neigh_e_prev2layer(index,hitEnergy[iN][jH]);
		if (newlayer == (newhitlayer+2)) myHit[iN].set_neigh_e_next2layer(index,hitEnergy[iN][jH]);
	      }
	    }
	  }//within radius
	}//within +/- 2 layers
	
      }//loop on hits to find neighbours

      for (unsigned iN(0); iN<nNoise;++iN){
	//select tracks
	if (myHit[iN].getMaxEnergy(-1)<threshBeforeAfterMin || myHit[iN].getMaxEnergy(1)<threshBeforeAfterMin || myHit[iN].getMaxEnergy(-1)>threshBeforeAfterMax || myHit[iN].getMaxEnergy(1)>threshBeforeAfterMax){
	  pass[iN]=false;	  
	}
	if (iN==0 && pass[0]) passSel[1]++;
	if (noise[iN]>0 && myHit[iN].getMaxEnergy(0)>threshCentralNeighMax*noise[iN]){
	  pass[iN]=false;	  
	} 
	if (pass[iN]) {
	  miphitvec[iN].push_back(myHit[iN]);
	  //std::cout << " vec pointer end loop " << iN << " " << &miphitvec[iN] << std::endl;
	}
	if (iN==0 && pass[0]) passSel[2]++;
      }//loop on noise
    }//loop on hits
    if (debug) std::cout << " --- Hits passing central Hit selection with 0 noise: " << passSel[0] << " +/-1 layers " << passSel[1] << " central neighbours " << passSel[2] << std::endl;

    for (unsigned iN(0); iN<nNoise;++iN){
      if (miphitvec[iN].size()>=nHitsInit) {
	std::cout << " -- problem with hits vector initial size, need bigger!"
		  << " Event " << ievt << " noise " << noise[iN] 
		  << " nHits = " << miphitvec[iN].size() << std::endl;
	//exit(0);
      }
    }

    miptree->Fill();

    geomConv.initialiseHistos();

  }//loop on entries


  outputFile->cd();
  miptree->Write();
  outputFile->Close();

  return 0;
  
}//main
