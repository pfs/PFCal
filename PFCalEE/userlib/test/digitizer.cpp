#include<string>
#include<set>
#include<iostream>
#include<fstream>
#include<sstream>
#include <boost/algorithm/string.hpp>
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "boost/function.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "TH2Poly.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "Math/Vector4D.h"

#include "TSystemDirectory.h"
#include "TSystemFile.h"

#include "fastjet/ClusterSequence.hh"

#include "HGCSSEvent.hh"
#include "HGCSSInfo.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSRecoJet.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDigitisation.hh"
#include "HGCSSDetector.hh"
#include "HGCSSGeometryConversion.hh"

using namespace fastjet;

using boost::lexical_cast;
namespace po=boost::program_options;

template <class T>
void extractParameterFromStr(std::string aStr,T & vec){ 
  if (aStr == "") return;
  std::vector<std::string> layVec;
  boost::split( layVec, aStr, boost::is_any_of(","));

  for (unsigned iE(0); iE<layVec.size(); ++iE){//loop on elements
    std::vector<std::string> lPair;
    boost::split( lPair, layVec[iE], boost::is_any_of(":"));
    if (lPair.size() != 2) {
      std::cout << " -- Wrong string for parameter given as input:" << layVec[iE] << " Try again, expecting exactly one symbol \":\" between two \",\" ..." << std::endl;
      exit(1);
    }
    std::vector<std::string> lLay;
    boost::split( lLay, lPair[0], boost::is_any_of("-"));
    if (lLay.size() > 2) {
      std::cout << " -- Wrong string for granularities given as input:" << lPair[0] << " Try again, expecting at most one symbol \"-\"." << std::endl;
      exit(1);
    }
    unsigned beginIdx =  atoi(lLay[0].c_str());
    unsigned endIdx = lLay.size() == 1 ? beginIdx :  atoi(lLay[1].c_str());
    for (unsigned iL(beginIdx); iL<endIdx+1; ++iL){
      if (iL < vec.size())
	std::istringstream(lPair[1])>>vec[iL];
      else {
	std::cout << " -- WARNING! Input parameter has more layers: " << endIdx << " than detector : " 
		  << vec.size()
		  << ". Ignoring additional layer #" << iL << "... PLEASE CHECK SETTINGS ARE CORRECT FOR EXISTING LAYERS!!"
		  << std::endl;
      }
    }
  }//loop on elements
}

/*
void processHist(const unsigned iL,
		 std::map<unsigned,MergeCells> & histE,
		 HGCSSDigitisation & myDigitiser,
		 TH1F* & p_noise,
		 const TH2Poly* histZ,
		 const double & meanZpos,
		 const bool isTBsetup,
		 const HGCSSSubDetector & subdet,
		 const std::vector<unsigned> & pThreshInADC,
		 const bool pSaveDigis,
		 HGCSSRecoHitVec & lDigiHits,
		 HGCSSRecoHitVec & lRecoHits,
		 const bool pMakeJets,
		 std::vector<PseudoJet> & lParticles
		 ){

  bool doSaturation=true;

  DetectorEnum adet = subdet.type;
  bool isScint = subdet.isScint;
  bool isSi = subdet.isSi;

  //double rLim = subdet.radiusLim;

  TIter next(histE->GetBins());
  TObject *obj=0; 
  TH2PolyBin *polyBin = 0;
  
  while ((obj=next())){
    polyBin=(TH2PolyBin*)obj;
    double x = (polyBin->GetXMax()+polyBin->GetXMin())/2.;
    double y = (polyBin->GetYMax()+polyBin->GetYMin())/2.;
    //double radius = sqrt(pow(x,2)+pow(y,2));
    
    double digiE = 0;
    double simE = polyBin->GetContent();
    //double time = 0;
    //if (simE>0) time = histTime->GetBinContent(iX,iY)/simE;
    
    //fill vector with neighbours and calculate cross-talk
    double xtalkE = simE;
    if (isScint){
      std::vector<double> simEvec;
      simEvec.push_back(simE);
      double side = polyBin->GetXMax()-polyBin->GetXMin();
      simEvec.push_back(histE->GetBinContent(histE->FindBin(x-side,y)));
      simEvec.push_back(histE->GetBinContent(histE->FindBin(x+side,y)));
      simEvec.push_back(histE->GetBinContent(histE->FindBin(x,y-side)));
      simEvec.push_back(histE->GetBinContent(histE->FindBin(x,y+side)));
      xtalkE = myDigitiser.ipXtalk(simEvec);
    }
      
    //bool passTime = myDigitiser.passTimeCut(adet,time);
    //if (!passTime) continue;
      
    double posz = meanZpos;
    //for noise only hits
    //if (simE>0) posz = histZ->GetBinContent(iX,iY)/simE;
    //else posz = meanZpos;
    
    
    //if (fabs(x) > 500 || fabs(y)>500) std::cout << " x=" << x << ", y=" << y << std::endl;
    
    
    //correct for particle angle in conversion to MIP
    //not necessary, if not done for aborber thickness either
    double simEcor = xtalkE;//isTBsetup ? xtalkE : myDigitiser.mipCor(xtalkE,x,y,posz);
    digiE = simEcor;
    
    if (isScint && simEcor>0 && doSaturation) {
      digiE = myDigitiser.digiE(simEcor);
    }
    myDigitiser.addNoise(digiE,iL,p_noise);
    
    double noiseFrac = 1.0;
    if (simEcor>0) noiseFrac = (digiE-simEcor)/simEcor;
    
    //for silicon-based Calo
    unsigned adc = 0;
    if (isSi){
      adc = myDigitiser.adcConverter(digiE,adet);
      digiE = myDigitiser.adcToMIP(adc,adet);
    }
    bool aboveThresh = //digiE > 0.5;
      (isSi && adc > pThreshInADC[iL]) ||
      (isScint && digiE > pThreshInADC[iL]*myDigitiser.adcToMIP(1,adet,false));
    //histE->SetBinContent(iX,iY,digiE);
    if ((!pSaveDigis && aboveThresh) ||
	pSaveDigis)
      {//save hits
	//double calibE = myDigitiser.MIPtoGeV(subdet,digiE);
	HGCSSRecoHit lRecHit;
	lRecHit.layer(iL);
	lRecHit.energy(digiE);
	lRecHit.adcCounts(adc);
	lRecHit.x(x);
	lRecHit.y(y);
	lRecHit.z(posz);
	lRecHit.noiseFraction(noiseFrac);
	//unsigned x_cell = static_cast<unsigned>(fabs(x)/(cellSize*granularity[iL]));
	//unsigned y_cell = static_cast<unsigned>(fabs(y)/(cellSize*granularity[iL]));
	//lRecHit.encodeCellId(x>0,y>0,x_cell,y_cell,granularity[iL]);
	
	if (pSaveDigis) lDigiHits.push_back(lRecHit);
	
	lRecoHits.push_back(lRecHit);
	
	if (pMakeJets){
	  if (posz>0) lParticles.push_back( PseudoJet(lRecHit.px(),lRecHit.py(),lRecHit.pz(),lRecHit.E()));
	}
	
      }//save hits
  }//loop on bins
     
}//processHist
*/
double getEtaPhys(const double & posx,
		  const double & posy,
		  const double & posz,
		  const HGCSSEvent * event){

    double posxcor = posx-event->vtx_x();
    double posycor = posy-event->vtx_y();
    double poszcor = posz-event->vtx_z();
    double r3d = sqrt(pow(posxcor,2)+pow(posycor,2)+pow(poszcor,2));
    double thetaPhys = acos(fabs(poszcor)/r3d);
    double etaPhys = (poszcor)>0? -1.*log(tan(thetaPhys/2.)) : log(tan(thetaPhys/2.));
    return etaPhys;
}

bool processEvent(const HGCSSEvent * event,
		  const std::vector<HGCSSGenParticle> * genvec,
		  const std::vector<HGCSSSimHit> * hitvec,
		  HGCSSGeometryConversion & geomConv,
		  HGCSSCalibration & mycalib,
		  HGCSSDetector & myDetector,
		  HGCSSDigitisation & myDigitiser,
		  const double etamean,
		  const double deta,
		  const bool doEtaSel,
		  const bool pSaveSims,
		  const unsigned debug,
		  const unsigned shape,
		  HGCSSGenParticleVec & lGenParts,
		  HGCSSSimHitVec & lSimHits
		  ){
  bool passGenFilter = false;

  if (doEtaSel){
    for (unsigned iP(0); iP<(*genvec).size(); ++iP){//loop on gen particles
      HGCSSGenParticle lPart = (*genvec)[iP];
      if (fabs(lPart.charge())>0.1 && fabs(lPart.eta()-etamean)<deta) {
	passGenFilter = true;
	lGenParts.push_back(lPart);
      }
    }
  }
  
  for (unsigned iH(0); iH<(*hitvec).size(); ++iH){//loop on hits
    HGCSSSimHit lHit = (*hitvec)[iH];
    if (lHit.energy()<=0) continue;
      
    //do not save hits with 0 energy...
    if (pSaveSims) lSimHits.push_back(lHit);
    
    unsigned layer = lHit.layer();
    const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
    DetectorEnum type = subdet.type;
    if (debug > 1) std::cout << " - layer " << layer << " " << subdet.name << " " << layer-subdet.layerIdMin << std::endl;


    std::pair<double,double> xy = lHit.get_xy(subdet,geomConv,shape);
    double posx = xy.first;//lHit.get_x(cellSize);
    double posy = xy.second;//lHit.get_y(cellSize);
    double radius = sqrt(pow(posx,2)+pow(posy,2));

    if (lHit.silayer() >= geomConv.getNumberOfSiLayers(type,radius)) continue;

    double posz = lHit.get_z();
    double etaPhys = getEtaPhys(posx,posy,posz,event);

    if (doEtaSel){
      bool passeta = fabs(etaPhys-etamean)<deta;
      if (!passeta) continue;
    }
    
    double energy = lHit.energy()*mycalib.MeVToMip(layer,radius); // if (energy > 0) std::cout << "sim energy = "<<lHit.energy()<<", reco energy = "<<energy<<std::endl;
    //mycalib already corrects for vertex position....
    double realtime = mycalib.correctTime(lHit.time(),posx,posy,posz);
    bool passTime = myDigitiser.passTimeCut(type,realtime);
    if (!passTime) continue;
    if (debug > 1) std::cout << " hit " << iH 
			     << " lay " << layer  
			     << " x " << posx 
			     << " y " << posy
			     << " z " << posz
			     << " t " << lHit.time() << " " << realtime
			     << std::endl;
    //geomConv.fill(type,subdetLayer,energy,realtime,posx,posy,posz);
    geomConv.fill(layer,energy,realtime,lHit.cellid(),posz);
    
  }//loop on input simhits

  return passGenFilter;
}


void processHist(const unsigned iL,
		 std::map<unsigned,MergeCells> & histE,
		 std::map<int,std::pair<double,double> > & geom,
		 HGCSSDigitisation & myDigitiser,
		 TH1F* & p_noise,
		 //const TH2Poly* histZ,
		 const double & meanZpos,
		 const bool isTBsetup,
		 const HGCSSSubDetector & subdet,
		 const std::vector<unsigned> & pThreshInADC,
		 const bool pSaveDigis,
		 HGCSSRecoHitVec & lDigiHits,
		 HGCSSRecoHitVec & lRecoHits,
		 const bool pMakeJets,
		 std::vector<PseudoJet> & lParticles
		 ){

  bool doSaturation=false;//true;

  DetectorEnum adet = subdet.type;
  bool isScint = subdet.isScint;
  bool isSi = subdet.isSi;
  //double rLim = subdet.radiusLim;
  std::map<unsigned,MergeCells>::iterator lIter = histE.begin();
  for (; lIter!=histE.end();++lIter){//loop on elements of the map
    //bin numbering starts at 1....
    //get bin number of map element iele
    unsigned iB = lIter->first;
    //cut overflows: very large IDs from the TH2Poly.
    if(iB>4000000000) continue;
    std::pair<double,double> xy = geom[iB];
    if (isScint) HGCSSGeometryConversion::convertFromEtaPhi(xy,meanZpos);
    double digiE = 0;
    double simE = lIter->second.energy;
    //double time = 0;
    //if (simE>0) time = histE[iele].time/simE;
    
    //fill vector with neighbours and calculate cross-talk
    double xtalkE = simE;
    //CAMM @TODO
    /*if (isScint){
      std::vector<double> simEvec;
      simEvec.push_back(simE);
      double side = polyBin->GetXMax()-polyBin->GetXMin();
      simEvec.push_back(histE->GetBinContent(histE->FindBin(x-side,y)));
      simEvec.push_back(histE->GetBinContent(histE->FindBin(x+side,y)));
      simEvec.push_back(histE->GetBinContent(histE->FindBin(x,y-side)));
      simEvec.push_back(histE->GetBinContent(histE->FindBin(x,y+side)));
      xtalkE = myDigitiser.ipXtalk(simEvec);
      }*/
      
    //bool passTime = myDigitiser.passTimeCut(adet,time);
    //if (!passTime) continue;
      
    double posz = meanZpos;
    //for noise only hits
    //if (simE>0) posz = histZ->GetBinContent(iX,iY)/simE;
    //else posz = meanZpos;
    
    
    //if (fabs(x) > 500 || fabs(y)>500) std::cout << " x=" << x << ", y=" << y << std::endl;
    
    
    //correct for particle angle in conversion to MIP
    //not necessary, if not done for aborber thickness either
    double simEcor = xtalkE;//isTBsetup ? xtalkE : myDigitiser.mipCor(xtalkE,x,y,posz);
    digiE = simEcor;
    
    if (isScint && simEcor>0 && doSaturation) {
      digiE = myDigitiser.digiE(simEcor);
    }
    myDigitiser.addNoise(digiE,iL,p_noise);
    
    double noiseFrac = 1.0;
    if (simEcor>0) noiseFrac = (digiE-simEcor)/simEcor;
    
    //for silicon-based Calo
    unsigned adc = 0;
    if (isSi){
      adc = myDigitiser.adcConverter(digiE,adet);
      digiE = myDigitiser.adcToMIP(adc,adet);
    }
    bool aboveThresh = //digiE > 0.5;
      (isSi && adc >= pThreshInADC[iL]) ||
      (isScint && digiE >= pThreshInADC[iL]*myDigitiser.adcToMIP(1,adet,false));
    //histE->SetBinContent(iX,iY,digiE);
    if ((!pSaveDigis && aboveThresh) ||
	pSaveDigis)
      {//save hits
	//double calibE = myDigitiser.MIPtoGeV(subdet,digiE);
	HGCSSRecoHit lRecHit;
	lRecHit.layer(iL);
	lRecHit.energy(digiE);
	lRecHit.adcCounts(adc);
	lRecHit.x(xy.first);
	lRecHit.y(xy.second);
	lRecHit.z(posz);
	lRecHit.noiseFraction(noiseFrac);
	//unsigned x_cell = static_cast<unsigned>(fabs(x)/(cellSize*granularity[iL]));
	//unsigned y_cell = static_cast<unsigned>(fabs(y)/(cellSize*granularity[iL]));
	//lRecHit.encodeCellId(x>0,y>0,x_cell,y_cell,granularity[iL]);
	
	if (pSaveDigis) lDigiHits.push_back(lRecHit);
	
	lRecoHits.push_back(lRecHit);
	
	if (pMakeJets){
	  if (posz>0) lParticles.push_back( PseudoJet(lRecHit.px(),lRecHit.py(),lRecHit.pz(),lRecHit.E()));
	}
	
      }//save hits
  }//loop on bins
     
}//processHist


int main(int argc, char** argv){//main  

  const unsigned evtmin = 0;//100;
  /////////////////////////////////////////////////////////////
  //parameters
  /////////////////////////////////////////////////////////////
  //Input output and config options
  std::string cfg;
  unsigned pNevts;
  std::string inFilePath;
  std::string outFilePath;
  std::string granulStr;//granularities layer_i-layer_j:factor,layer:factor,...
  std::string noiseStr;//noise (in Mips) layer_i-layer_j:factor,layer:factor,...
  std::string threshStr;//threshold (in ADC counts) layer_i-layer_j:factor,layer:factor,...
  unsigned interCalib;//intercalib factor in %
  unsigned nSiLayers;//Number of si layers for TB setups
  unsigned nPU;//number of PU to overlay
  std::string puPath;
  //for selecting a ring in eta - for noise studies.
  double etamean;//ring etamean (default=0=no eta sel)
  double deta;//ring +/- delta eta value
  //remove noise hits everywhere to speed up processing, e.g. for muons.
  bool addNoiseHits;

  unsigned pSeed ;
  unsigned debug;
  bool pSaveDigis;
  bool pSaveSims;
  bool pMakeJets;

  bool pFilterOnGenParticles;

  po::options_description preconfig("Configuration"); 
  preconfig.add_options()("cfg,c",po::value<std::string>(&cfg)->required());
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(preconfig).allow_unregistered().run(), vm);
  po::notify(vm);
  po::options_description config("Configuration");
  config.add_options()
    ("pNevts,n",      po::value<unsigned>(&pNevts)->default_value(0))
    ("inFilePath,i",  po::value<std::string>(&inFilePath)->required())
    ("outFilePath,o", po::value<std::string>(&outFilePath)->required())
    ("granulStr",     po::value<std::string>(&granulStr)->required())
    ("noiseStr",      po::value<std::string>(&noiseStr)->required())
    ("threshStr",     po::value<std::string>(&threshStr)->required())
    ("interCalib",    po::value<unsigned>(&interCalib)->default_value(3))
    ("nSiLayers",     po::value<unsigned>(&nSiLayers)->default_value(2))
    ("nPU",           po::value<unsigned>(&nPU)->default_value(0))
    ("puPath",        po::value<std::string>(&puPath)->default_value(""))
    ("etamean",       po::value<double>(&etamean)->default_value(0))
    ("deta",          po::value<double>(&deta)->default_value(10))
    ("addNoiseHits,a",po::value<bool>(&addNoiseHits)->default_value(true))
    ("pSeed,s",       po::value<unsigned>(&pSeed)->default_value(0))
    ("debug,d",       po::value<unsigned>(&debug)->default_value(0))
    ("pSaveDigis",    po::value<bool>(&pSaveDigis)->default_value(false))
    ("pSaveSims",     po::value<bool>(&pSaveSims)->default_value(false))
    ("pMakeJets",     po::value<bool>(&pMakeJets)->default_value(false))
    ("pFilterOnGenParticles",     po::value<bool>(&pFilterOnGenParticles)->default_value(true))
    ;

  po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
  po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
  po::notify(vm);




  if (nPU>0 && puPath.size()==0) {
    std::cout << " -- Error! Missing full path to minbias file. Exiting." << std::endl;
    return 1;
  }

  std::cout << " ----------------------------------------" << std::endl
            << " -- Input parameters: " << std::endl
            << " -- Input file path: " << inFilePath << std::endl
            << " -- Output file path: " << outFilePath << std::endl
            << " -- Processing " ;
  if (pNevts>0) std::cout << pNevts;
  else std::cout << "all";
  std::cout << " events." << std::endl
            << " -- Granularities: " << granulStr << std::endl
            << " -- noise: " << noiseStr << std::endl
            << " -- thresholds: " << threshStr << std::endl
            << " -- intercalibration factor (in %): " << interCalib << std::endl
	    << " -- number of Si layers: " << nSiLayers << std::endl
            << " -- number of PU: " << nPU << std::endl
	    << " -- pu file path: " << puPath << std::endl
    ;


  bool doEtaSel = false;
  if (etamean > 1.4){
    std::cout << " -- Eta selection: " << etamean << " +/- " << deta << std::endl;
    doEtaSel = true;
  }

  if (debug>0) std::cout << " -- DEBUG output is set to " << debug << std::endl;
  
  //try to get model automatically
  //if (inFilePath.find("model0")!=inFilePath.npos) pModel = "model0";
  //else if (inFilePath.find("model1")!=inFilePath.npos) pModel = "model1";
  //else if (inFilePath.find("model2")!=inFilePath.npos) pModel = "model2";
  //else if (inFilePath.find("model3")!=inFilePath.npos) pModel = "model3";


  //std::cout << " -- Model is set to : " << pModel << std::endl;
  std::cout<< " -- Random seed will be set to : " << pSeed << std::endl;
  if (pSaveDigis) std::cout << " -- DigiHits are saved." << std::endl;
  if (pSaveSims) std::cout << " -- SimHits are saved." << std::endl;
  if (pMakeJets) std::cout << " -- Making jets." << std::endl;
  std::cout << " ----------------------------------------" << std::endl;
  
  //////////////////////////////////////////////////////////
  //// Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////
  //for HGCAL, true means only 12 FHCAL layers considered (24 are simulated)
  bool concept = true;

  // choose a jet definition
  double R = 0.5;
  JetDefinition jet_def(antikt_algorithm, R);

  // define the outer edge of the scint layers
  double outerScintBoundary [69] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				     0, 0, 0,
				     1.474, 1.455, 1.437, 1.420, 1.403, 1.376, 1.351, 1.327, 
				     1.306, 1.316, 1.332, 1.348, 1.364, 1.379, 1.395, 1.410};



  //////////////////////////////////////////////////////////
  //// End Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////
  //input
  /////////////////////////////////////////////////////////////

  std::string inputStr = inFilePath ;
  TFile *inputFile = TFile::Open(inputStr.c_str());

  std::cout << "Opening input file " << inputStr << std::flush;

  if (!inputFile) {
    std::cout << " -- Error, input file " << inputStr << " cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  std::cout << "...Done." << std::endl <<  "Getting HGCSSTree" << std::flush;

  TTree *inputTree = (TTree*)inputFile->Get("HGCSSTree");
  if (!inputTree){
    std::cout << " -- Error, tree HGCSSTree  cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  std::cout << "...Done." <<  std::endl;
  /////////////////////////////////////////////////////////////
  //input tree
  /////////////////////////////////////////////////////////////

  HGCSSInfo * info=(HGCSSInfo*)inputFile->Get("Info");

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
  if (isTBsetup) std::cout << " -- Number of Si layers: " << nSiLayers << std::endl;
  else std::cout << " -- Number of Si layers ignored: hardcoded as a function of radius in HGCSSGeometryConversion class." << std::endl;


  /////////////////////////////////////////////////////////////
  //  //input PU hits
  /////////////////////////////////////////////////////////////

  bool signalIsPu = false;

  TChain *puTree = new TChain("HGCSSTree");
  unsigned nPuVtx = 0;
  unsigned nPuEvts = 0;
  std::vector<HGCSSSimHit> * puhitvec = 0;
  std::vector<HGCSSGenParticle> * pugenvec = 0;

  if(nPU!=0){
    
    TString localMountPuPath(puPath.c_str());
    localMountPuPath.ReplaceAll("root://eoscms/","");
    if(localMountPuPath.BeginsWith("/store")) localMountPuPath="/eos/cms"+localMountPuPath;
    
    TSystemDirectory dir(localMountPuPath.Data(),localMountPuPath.Data());
    TList *files = dir.GetListOfFiles();
    TSystemFile *file;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      if( file->IsDirectory() ) continue;
      TString fname( file->GetName() );
      TString lversion = "version";
      lversion += versionNumber;
      if( !fname.Contains(".root") ) continue;
      if( !fname.Contains("HGcal") ) continue;
      if( !fname.Contains(lversion) ) continue;
      if(  fname.Contains("Digi") ) continue;
      TString puInput(localMountPuPath+"/"+fname);
      if (puInput==inFilePath.c_str()) {
        std::cout << " -- duplicate: file already used as signal. Removing." << std::endl; 
        signalIsPu = true;
        continue;
      }
      puTree->AddFile(puInput);
      std::cout << "Adding MinBias file:" << puInput << std::endl;
    }

    puTree->SetBranchAddress("HGCSSSimHitVec",&puhitvec);
    puTree->SetBranchAddress("HGCSSGenParticleVec",&pugenvec);
    nPuEvts = puTree->GetEntries();
    std::cout << "- Number of PU events available: " << nPuEvts  << std::endl;
  }
  /////////////////////////////////////////////////////////////
  //input signal tree
  /////////////////////////////////////////////////////////////

  HGCSSEvent * event=0;
  std::vector<HGCSSSimHit> * hitvec = 0;
  std::vector<HGCSSGenParticle> * genvec = 0;

  inputTree->SetBranchAddress("HGCSSEvent",&event);
  inputTree->SetBranchAddress("HGCSSSimHitVec",&hitvec);
  inputTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);

  //initialise detector
  HGCSSDetector & myDetector = theDetector();
  bool isCaliceHcal = versionNumber==23;//inFilePath.find("version23")!=inFilePath.npos || inFilePath.find("version_23")!=inFilePath.npos;

  //unsigned versionNumber = 0;
  //unsigned nchar = 0;
  //if (inFilePath.substr(inFilePath.find("version")+8,1)=="_") nchar = 1;
  //else nchar = 2;
  //std::string lvers = inFilePath.substr(inFilePath.find("version")+7,nchar);
  //std::istringstream(lvers)>>versionNumber;
  
  std::cout << " -- Calor size XY = " << calorSizeXY
	    << ", version number = " << versionNumber 
	    << ", model = " << model
	    << ", cellSize = " << cellSize
	    << ", shape = " << shape
	    << std::endl;

  bool bypassR = false;
  if (isTBsetup) bypassR = true;
  myDetector.buildDetector(versionNumber,model,concept,isCaliceHcal,bypassR);

  //initialise calibration class
  HGCSSCalibration mycalib(inFilePath,bypassR,nSiLayers);

  const unsigned nLayers = myDetector.nLayers();

  HGCSSGeometryConversion geomConv(model,cellSize,bypassR,nSiLayers);
  geomConv.setXYwidth(calorSizeXY);
  geomConv.setVersion(versionNumber);
  //const double xWidth = geomConv.getXYwidth();

  if (shape==2) geomConv.initialiseDiamondMap(calorSizeXY,10.);
  else if (shape==3) geomConv.initialiseTriangleMap(calorSizeXY,10.*sqrt(2.));
  else if (shape==1) geomConv.initialiseHoneyComb(calorSizeXY,cellSize);
  else if (shape==4) geomConv.initialiseSquareMap(calorSizeXY,10.);

  //square map for BHCAL
  //geomConv.initialiseSquareMap1(1.4,3.0,-1.*TMath::Pi(),TMath::Pi(),0.01745);//eta phi segmentation
  //geomConv.initialiseSquareMap2(1.4,3.0,-1.*TMath::Pi(),TMath::Pi(),0.02182);//eta phi segmentation
  geomConv.initialiseSquareMap1(1.3,3.0,-1.*TMath::Pi(),TMath::Pi(),TMath::Pi()*2./360.);//eta phi segmentation, hardcoded '1.3' to include outer edge fo BH
  geomConv.initialiseSquareMap2(1.3,3.0,-1.*TMath::Pi(),TMath::Pi(),TMath::Pi()*2./288.);//eta phi segmentation, hardcoded '1.3' to include outer edge fo BH
  //geomConv.initialiseSquareMap2(1.4,3.0,0,2*TMath::Pi(),0.02618);//eta phi segmentation

  HGCSSDigitisation myDigitiser;
  myDigitiser.setIntercalibrationFactor(interCalib);
  std::cout << " -- Intercalibration factor set to (%): " << interCalib << std::endl;

  std::vector<unsigned> granularity;
  granularity.resize(nLayers,1);
  std::vector<double> pNoiseInMips;
  pNoiseInMips.resize(nLayers,0.12);
  std::vector<unsigned> pThreshInADC;
  pThreshInADC.resize(nLayers,5);

  extractParameterFromStr<std::vector<unsigned> >(granulStr,granularity);
  extractParameterFromStr<std::vector<double> >(noiseStr,pNoiseInMips);
  extractParameterFromStr<std::vector<unsigned> >(threshStr,pThreshInADC);

  //unsigned nbCells = 0;

  if (doEtaSel){
    for (unsigned iL(0); iL<nLayers; ++iL){
      pNoiseInMips[iL] = 0;
      pThreshInADC[iL] = 0;
    }
  }

  std::cout << " -- Granularities and noise are setup like this:" << std::endl;
  for (unsigned iL(0); iL<nLayers; ++iL){
    std::cout << "Layer " ;
    if (iL<10) std::cout << " ";
    std::cout << iL << " : " << granularity[iL] << ", " << pNoiseInMips[iL] << " mips, " << pThreshInADC[iL] << " adc - ";
    if (iL%5==4) std::cout << std::endl;
    myDigitiser.setNoise(iL,pNoiseInMips[iL]);
    //nbCells += N_CELLS_XY_MAX/(granularity[iL]*granularity[iL]);
  }
  std::cout << std::endl;     
  //std::cout << " -- Total number of cells = " << nbCells << std::endl;

  geomConv.setGranularity(granularity);
  geomConv.initialiseHistos();

  TRandom3 *lRndm = new TRandom3();
  lRndm->SetSeed(pSeed);
  myDigitiser.setRandomSeed(pSeed);

  std::cout << " -- Random3 seed = " << lRndm->GetSeed() << std::endl
	    << " ----------------------------------------" << std::endl;


  /////////////////////////////////////////////////////////////
  //output
  /////////////////////////////////////////////////////////////

  std::ostringstream outputStr;
  outputStr << outFilePath << "/DigiPFcal" ;
  //if (doEtaSel) {
  //outputStr << "_eta" << (etamean-deta) << "_" << (etamean+deta);
  //}
  if (pSaveDigis)  outputStr << "_withDigiHits";
  if (pSaveSims)  outputStr << "_withSimHits";
  outputStr << ".root";
  
  TFile *outputFile = TFile::Open(outputStr.str().c_str(),"RECREATE");

  if (!outputFile) {
    std::cout << " -- Error, output file " << outputStr.str() << " cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- File will be saved as " << outputStr.str() << std::endl;
  }

  HGCSSInfo *lInfo = new HGCSSInfo();
  lInfo->calorSizeXY(calorSizeXY);
  lInfo->cellSize(cellSize);
  lInfo->version(versionNumber);
  lInfo->model(model);
  lInfo->shape(shape);
  TTree *outputTree = new TTree("RecoTree","HGC Standalone simulation reco tree");
  HGCSSSimHitVec lSimHits;
  HGCSSRecoHitVec lDigiHits;
  HGCSSRecoHitVec lRecoHits;
  HGCSSRecoJetVec lCaloJets;
  HGCSSGenParticleVec lGenParts;
  unsigned maxSimHits = 0;
  unsigned maxRecHits = 0;
  unsigned maxRecJets = 0;
  unsigned maxGenParts = 0;
  HGCSSEvent lEvent;
  outputTree->Branch("HGCSSEvent",&lEvent);
  if (nPU!=0) outputTree->Branch("nPuVtx",&nPuVtx);
  if (pSaveSims) outputTree->Branch("HGCSSSimHitVec","std::vector<HGCSSSimHit>",&lSimHits);
  if (pSaveDigis) outputTree->Branch("HGCSSDigiHitVec","std::vector<HGCSSRecoHit>",&lDigiHits);
  outputTree->Branch("HGCSSRecoHitVec","std::vector<HGCSSRecoHit>",&lRecoHits);
  outputTree->Branch("HGCSSGenParticleVec","std::vector<HGCSSGenParticle>",&lGenParts);

  if (pMakeJets) outputTree->Branch("HGCSSRecoJetVec","std::vector<HGCSSRecoJet>",&lCaloJets);
  TH1F * p_noise = new TH1F("noiseCheck",";noise (MIPs)",100,-5,5);


  /////////////////////////////////////////////////////////////
  //Loop on events
  /////////////////////////////////////////////////////////////

  const unsigned nEvts = (pNevts > inputTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(inputTree->GetEntries()) : pNevts;

  std::cout << "- Processing = " << nEvts  << " events out of " << inputTree->GetEntries() << std::endl;

  std::vector<PseudoJet> lParticles;

  for (unsigned ievt(evtmin); ievt<evtmin+nEvts; ++ievt){//loop on entries
    inputTree->GetEntry(ievt);

    lEvent.eventNumber(event->eventNumber());
    lEvent.vtx_x(event->vtx_x());
    lEvent.vtx_y(event->vtx_y());
    lEvent.vtx_z(event->vtx_z());
    //unsigned layer = volNb;
    
    mycalib.setVertex(lEvent.vtx_x(),lEvent.vtx_y(),lEvent.vtx_z());


    if (debug>0) {
      std::cout << " **DEBUG** Processing evt " << ievt << std::endl;
    }
    else if (ievt%50 == 0) std::cout << "... Processing event: " << ievt << std::endl;
    
    bool passGenFilter = processEvent(event,genvec,hitvec,geomConv,mycalib,myDetector,myDigitiser,etamean,deta,doEtaSel,pSaveSims,debug,shape,lGenParts,lSimHits);

    if(nPU!=0){
      //get PU events
      //std::vector<unsigned> ipuevt;

      //get poisson <140>
      nPuVtx = lRndm->Poisson(nPU);
      //ipuevt.resize(nPuVtx,1);
      if (signalIsPu) nPuVtx -= 1;

      //pick first event random then sequencial.... better for disk access
      unsigned ipuevt = 0;
      while (1){
	ipuevt = lRndm->Integer(nPuEvts);
	if ( (ipuevt+nPuVtx) < nPuEvts) break;
      }

      std::cout << " -- Adding " << nPuVtx << " events to signal event: " << ievt << " starting from index " << ipuevt << std::endl;

      //change event by event: bit problematic for disk access....
      /*
	std::set<unsigned> lidxSet;
	for (unsigned iV(0); iV<nPuVtx; ++iV){//loop on interactions
        unsigned ipuevt = 0;
	while (1){
	ipuevt = lRndm->Integer(nPuEvts);
	if (lidxSet.find(ipuevt)==lidxSet.end()){
	lidxSet.insert(ipuevt);
	break;
	}
	else {
	std::cout << " -- Found duplicate ! Taking another shot." << std::endl;
	}
	}
	//std::cout << " ---- adding evt " << ipuevt << std::endl;
	*/
      for (unsigned iV(ipuevt); iV<ipuevt+nPuVtx; ++iV){//loop on interactions
        puTree->GetEntry(iV);
	
	passGenFilter = passGenFilter || 
	  processEvent(event,pugenvec,puhitvec,geomConv,mycalib,myDetector,myDigitiser,etamean,deta,doEtaSel,pSaveSims,debug,shape,lGenParts,lSimHits);
		
      }//loop on interactions
    }//add PU

    if (debug>0) {
      std::cout << " **DEBUG** simhits = " << (*hitvec).size() << " " << lSimHits.size() << std::endl;
    }

    if ((doEtaSel && 
	 ((pFilterOnGenParticles && passGenFilter)||!pFilterOnGenParticles)) ||
	!doEtaSel) {//filter on genParts if etasel

      //create hits, everywhere to have also pure noise
      //digitise
      //apply threshold
      //save
      unsigned nTotBins = 0;
      for (unsigned iL(0); iL<nLayers; ++iL){//loop on layers
	//TH2Poly *histE = geomConv.get2DHist(iL,"E");
	std::map<unsigned,MergeCells> & histE = geomConv.get2DHist(iL);
	//TH2Poly *histTime = geomConv.get2DHist(iL,"Time");
	//TH2Poly *histZ = geomConv.get2DHist(iL,"Z");
	const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(iL);
	bool isScint = subdet.isScint;
	//nTotBins += histE->GetNumberOfBins();
	
	std::map<int,std::pair<double,double> > & geom = isScint?(subdet.type==DetectorEnum::BHCAL1?geomConv.squareGeom1:geomConv.squareGeom2): shape==4?geomConv.squareGeom:shape==2?geomConv.diamGeom:shape==3?geomConv.triangleGeom:geomConv.hexaGeom;
	
	unsigned nBins = geom.size();//isScint||shape==4?geomConv.squareMap()->GetNumberOfBins() : shape==2?geomConv.diamondMap()->GetNumberOfBins() : shape==3? geomConv.triangleMap()->GetNumberOfBins() : geomConv.hexagonMap()->GetNumberOfBins();
	nTotBins += nBins;
	if (pSaveDigis) lDigiHits.reserve(nTotBins);
	
	//double meanZpos = geomConv.getAverageZ(iL);
	double meanZpos = myDetector.sensitiveZ(iL);
	double etaBoundary = myDetector.etaBoundary(iL);
	//extend map to include all cells in eta=1.4-3 region
	//in eta ring if saving only one eta ring....
	if (addNoiseHits) {
	  for (unsigned iB(1); iB<nBins+1;++iB){
	    std::pair<double,double> xy = geom[iB];
	    if (isScint) {
	      HGCSSGeometryConversion::convertFromEtaPhi(xy,meanZpos);
	    }

	    double etaPhys = getEtaPhys(xy.first,xy.second,meanZpos,event);

	    ROOT::Math::XYZPoint lpos = ROOT::Math::XYZPoint(xy.first,xy.second,meanZpos);
	    double etaDet = lpos.eta();
	    bool passeta = etaDet>1.3 && etaDet<3.0;
	    if (doEtaSel) passeta = fabs(etaPhys-etamean)<deta;
	    else {
	      if (isScint) passeta = etaDet>outerScintBoundary[iL] && etaDet<=etaBoundary; // only simulate noise within the physical bounds of the detector
	      else passeta = etaDet>etaBoundary && etaDet<3.0;
	    }
	    if (!passeta) continue;
	    MergeCells tmpCell;
	    tmpCell.energy = 0;
	    tmpCell.time = 0;
	    histE.insert(std::pair<unsigned,MergeCells>(iB,tmpCell));
	  }
	}
	
	//std::cout << iL << " " << meanZpos << " map size " << histE.size() << std::endl;
	
	if (debug>0){
	  std::cout << " -- Layer " << iL << " " << subdet.name << " z=" << meanZpos
		    << " bins = " << nBins << " histE entries = " << histE.size() << std::endl;
	}
	
	//cell-to-cell cross-talk for scintillator
	if (isScint){
	  //2.5% per 30-mm edge
	  //myDigitiser.setIPCrossTalk(0.025*geomConv.cellSize(iL,0)/30.);
	}
	else {
	  myDigitiser.setIPCrossTalk(0);
	}
	
	//processHist(iL,histE,myDigitiser,p_noise,histZ,meanZpos,isTBsetup,subdet,pThreshInADC,pSaveDigis,lDigiHits,lRecoHits,pMakeJets,lParticles);
	
	processHist(iL,histE,geom,myDigitiser,p_noise,meanZpos,isTBsetup,subdet,pThreshInADC,pSaveDigis,lDigiHits,lRecoHits,pMakeJets,lParticles);
	
      }//loop on layers
      
      if (debug) {
	std::cout << " **DEBUG** sim-digi-reco hits = " << (*hitvec).size() << "-" << lDigiHits.size() << "-" << lRecoHits.size() << std::endl;
      }
      
      
      if (pMakeJets){//pMakeJets
	
	// run the clustering, extract the jets
	ClusterSequence cs(lParticles, jet_def);
	std::vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
	
	// print the jets
	std::cout <<   "-- evt " << ievt << ": found " << jets.size() << " Jets." << std::endl;
	for (unsigned i = 0; i < jets.size(); i++) {
	  const PseudoJet & lFastJet = jets[i];
	  //TOFIX // inverted y and z...
	  HGCSSRecoJet ljet(lFastJet.px(),
			    lFastJet.py(),
			    lFastJet.pz(),
			    lFastJet.E());
	  if (lFastJet.has_constituents()) ljet.nConstituents(lFastJet.constituents().size());
	  if (lFastJet.has_area()){
	    ljet.area(lFastJet.area());
	    ljet.area_error(lFastJet.area_error());
	  }
	  
	  lCaloJets.push_back(ljet);
	  std::cout << " -------- jet " << i << ": "
		    << lFastJet.E() << " " 
		    << lFastJet.perp() << " " 
		    << lFastJet.rap() << " " << lFastJet.phi() << " "
		    << lFastJet.constituents().size() << std::endl;
	  // std::vector<PseudoJet> constituents = lFastJet.constituents();
	  // for (unsigned j = 0; j < constituents.size(); j++) {
	  //   std::cout << "    constituent " << j << "'s pt: " << constituents[j].perp()
	  // 	      << std::endl;
	  // }
	}
	
      }//pMakeJets
      
      outputTree->Fill();

    }//filter on genParts if etasel

    //reserve necessary space and clear vectors.
    if (lSimHits.size() > maxSimHits) {
      maxSimHits = 2*lSimHits.size();
      std::cout << " -- INFO: event " << ievt << " maxSimHits updated to " << maxSimHits << std::endl;
    }
    if (lRecoHits.size() > maxRecHits) {
      maxRecHits = 2*lRecoHits.size();
      std::cout << " -- INFO: event " << ievt << " maxRecHits updated to " << maxRecHits << std::endl;
    }
    if (lCaloJets.size() > maxRecJets) {
      maxRecJets = 2*lCaloJets.size();
      std::cout << " -- INFO: event " << ievt << " maxRecJets updated to " << maxRecJets << std::endl;
    }
    if (lGenParts.size() > maxGenParts) {
      maxGenParts = 2*lGenParts.size();
      std::cout << " -- INFO: event " << ievt << " maxGenParts updated to " << maxGenParts << std::endl;
    }
    lSimHits.clear();
    lDigiHits.clear();
    lRecoHits.clear();
    lCaloJets.clear();
    lGenParts.clear();
    geomConv.initialiseHistos();
    lParticles.clear();
    if (pSaveSims) lSimHits.reserve(maxSimHits);
    lRecoHits.reserve(maxRecHits);
    lParticles.reserve(maxRecHits);
    lCaloJets.reserve(maxRecJets);
    lGenParts.reserve(maxGenParts);
  }//loop on entries

  outputFile->cd();
  outputFile->WriteObjectAny(lInfo,"HGCSSInfo","Info");
  outputTree->Write();
  p_noise->Write();
  outputFile->Close();

  return 0;

}//main
