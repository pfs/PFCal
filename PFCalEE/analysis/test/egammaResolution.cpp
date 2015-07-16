#include<string>
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

#include "PositionFit.hh"
#include "SignalRegion.hh"

#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"

using boost::lexical_cast;
namespace po=boost::program_options;

bool testInputFile(std::string input, TFile* & file){
  file = TFile::Open(input.c_str());
  
  if (!file) {
    std::cout << " -- Error, input file " << input.c_str() << " cannot be opened. Skipping..." << std::endl;
    return false;
  }
  else std::cout << " -- input file " << file->GetName() << " successfully opened." << std::endl;
  return true;
};


int main(int argc, char** argv){//main  

  //Input output and config options
  std::string cfg;
  bool concept;
  //size of signal region to perform Chi2 position fit.
  //in units of 2.5mm cells to accomodate different granularities
  unsigned nSR;
  //maximum value of residuals to use in error matrix: discard positions that are too far away 
  double residualMax;//mm
  unsigned pNevts;
  std::string filePath;
  std::string digifilePath;
  unsigned nRuns;
  std::string simFileName;
  std::string recoFileName;
  std::string outPath;
  unsigned nSiLayers;
  //0:do just the energies, 1:do fit+energies, 2: do zpos+fit+energies
  unsigned redoStep;
  unsigned debug;
  bool applyPuMixFix;

  po::options_description preconfig("Configuration"); 
  preconfig.add_options()("cfg,c",po::value<std::string>(&cfg)->required());
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(preconfig).allow_unregistered().run(), vm);
  po::notify(vm);
  po::options_description config("Configuration");
  config.add_options()
    //Input output and config options //->required()
    ("concept",        po::value<bool>(&concept)->default_value(true))
    ("nSR",            po::value<unsigned>(&nSR)->default_value(3))
    ("residualMax",    po::value<double>(&residualMax)->default_value(25))
    ("pNevts,n",       po::value<unsigned>(&pNevts)->default_value(0))
    ("filePath,i",     po::value<std::string>(&filePath)->required())
    ("digifilePath", po::value<std::string>(&digifilePath)->default_value(""))
    ("nRuns",        po::value<unsigned>(&nRuns)->default_value(0))
    ("simFileName,s",  po::value<std::string>(&simFileName)->required())
    ("recoFileName,r", po::value<std::string>(&recoFileName)->required())
    ("outPath,o",      po::value<std::string>(&outPath)->required())
    ("nSiLayers",      po::value<unsigned>(&nSiLayers)->default_value(2))
    ("redoStep",       po::value<unsigned>(&redoStep)->default_value(0))
    ("debug,d",        po::value<unsigned>(&debug)->default_value(0))
    ("applyPuMixFix",  po::value<bool>(&applyPuMixFix)->default_value(false))
    ;

  // ("output_name,o",            po::value<std::string>(&outputname)->default_value("tmp.root"))

  po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
  po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
  po::notify(vm);



  std::string inFilePath = filePath+simFileName;

  size_t end=outPath.find_last_of(".");
  std::string outFolder = outPath.substr(0,end);


  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << filePath << std::endl
	    << " -- Digi Input file path: " << digifilePath << std::endl
	    << " -- Output file path: " << outPath << std::endl
	    << " -- Output folder: " << outFolder << std::endl
	    << " -- Requiring " << nSiLayers << " si layers." << std::endl
	    << " -- Number cells in signal region for fit: " << nSR << " cells" << std::endl
	    << " -- Residual max considered for filling matrix and fitting: " << residualMax << " mm" << std::endl
	    << " -- Apply PUMix fix? " << applyPuMixFix << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events." << std::endl;

  TRandom3 lRndm(1);
  std::cout << " -- Random number seed: " << lRndm.GetSeed() << std::endl;

  /////////////////////////////////////////////////////////////
  //input
  /////////////////////////////////////////////////////////////

  std::ostringstream inputsim;
  inputsim << filePath << "/" << simFileName;
  std::ostringstream inputrec;
  if (digifilePath.size()==0)
    inputrec << filePath << "/" << recoFileName;
  else 
    inputrec << digifilePath << "/" << recoFileName;

  //std::cout << inputsim.str() << " " << inputrec.str() << std::endl;

  HGCSSInfo * info;

  TChain *lSimTree = new TChain("HGCSSTree");
  TChain *lRecTree = 0;
  
  TFile * simFile = 0;
  TFile * recFile = 0;

  if (recoFileName.find("Digi") != recoFileName.npos) 
    lRecTree = new TChain("RecoTree");
  else lRecTree = new TChain("PUTree");

  if (nRuns == 0){
    if (!testInputFile(inputsim.str(),simFile)) return 1;
    lSimTree->AddFile(inputsim.str().c_str());
    if (!testInputFile(inputrec.str(),recFile)) return 1;
    lRecTree->AddFile(inputrec.str().c_str());
  }
  else {
    for (unsigned i(0);i<nRuns;++i){
      std::ostringstream lstr;
      lstr << inputsim.str() << "_run" << i << ".root";
      if (testInputFile(lstr.str(),simFile)){  
	if (simFile) info =(HGCSSInfo*)simFile->Get("Info");
	else {
	  std::cout << " -- Error in getting information from simfile!" << std::endl;
	  return 1;
	}
      }
      else continue;
      std::ostringstream lstrrec;
      lstrrec << inputrec.str() << "_run" << i << ".root";
      if (!testInputFile(lstrrec.str(),recFile)) continue;
      //add only if both sim and reco are found...
      lSimTree->AddFile(lstr.str().c_str());
      lRecTree->AddFile(lstrrec.str().c_str());
    }
  }

  if (!lSimTree){
    std::cout << " -- Error, tree HGCSSTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  if (!lRecTree){
    std::cout << " -- Error, tree RecoTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }


  /////////////////////////////////////////////////////////////
  //Info
  /////////////////////////////////////////////////////////////

  const double cellSize = info->cellSize();
  const unsigned versionNumber = info->version();
  const unsigned model = info->model();
  
  //models 0,1 or 3.
  //bool isTBsetup = (model != 2);
  bool isCaliceHcal = versionNumber==23;//inFilePath.find("version23")!=inFilePath.npos || inFilePath.find("version_23")!=inFilePath.npos;

  //extract input energy

  std::cout << " -- Version number is : " << versionNumber 
	    << ", model = " << model
	    << ", cellSize = " << cellSize
	    << std::endl;


  //initialise detector
  HGCSSDetector & myDetector = theDetector();
 
  myDetector.buildDetector(versionNumber,concept,isCaliceHcal);

  const unsigned nLayers = myDetector.nLayers();
  const unsigned nSections = myDetector.nSections();

  std::cout << " -- N layers = " << nLayers << std::endl
	    << " -- N sections = " << nSections << std::endl;


  HGCSSGeometryConversion geomConv(inFilePath,model,cellSize);
  //set granularity to get cellsize for PU subtraction
  std::vector<unsigned> granularity;
  granularity.resize(nLayers,4);
  geomConv.setGranularity(granularity);

  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  ///////// Output File // /////////////////////////
  //////////////////////////////////////////////////
  //////////////////////////////////////////////////

  TFile *outputFile = TFile::Open(outPath.c_str(),"RECREATE");
  
  if (!outputFile) {
    std::cout << " -- Error, output file " << outPath << " cannot be opened. Please create output directory. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
  }
  outputFile->cd();


  ///initialise PU density object

  HGCSSPUenergy puDensity("data/EnergyDensity.dat");


    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
    ///////// positionFit /////////////////////////////
    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
  
  const unsigned nEvts = ((pNevts > lSimTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lSimTree->GetEntries()) : pNevts) ;
  

  bool doMatrix = false;

  PositionFit lChi2Fit(nSR,residualMax,nLayers,nSiLayers,applyPuMixFix,debug,doMatrix);
  lChi2Fit.initialise(outputFile,"PositionFit",outFolder,geomConv,puDensity);

  //try getting z position from input file, if doesn't exit,
  //perform first loop over simhits to find z positions of layers
  if ((redoStep<2 && !lChi2Fit.getZpositions(versionNumber)) || redoStep>1)
    lChi2Fit.getZpositions(versionNumber,lSimTree,nEvts);
  
  //perform second loop over events to find positions to fit and get energies
  SignalRegion SignalEnergy(outFolder, nLayers, nEvts, geomConv, puDensity,applyPuMixFix,versionNumber);
  SignalEnergy.initialise(outputFile,"Energies");

  //initialise
  bool dofit = redoStep>0 || !SignalEnergy.initialiseFitPositions();
  if (!dofit && redoStep==0) std::cout << " -- Info: fit positions taken from file on disk." << std::endl;
  else std::cout << " -- Info: redoing least square fit." << std::endl;

  if (redoStep>0 || (redoStep==0 && dofit)){
    if (doMatrix){
      lChi2Fit.getInitialPositions(lSimTree,lRecTree,nEvts);
      lChi2Fit.finaliseErrorMatrix();
    } else {
      lChi2Fit.initialiseClusterHistograms();
      lChi2Fit.initialisePositionHistograms();
      //lChi2Fit.setMatrixFolder(outFolder);
      lChi2Fit.setMatrixFolder("/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalDescopPosCheckBack/gitV02-01-01/version12/gamma/v5_30_/eta16_et70_pu0_IC3/");
      if (!lChi2Fit.fillMatrixFromFile(false)) {
	std::cout << " Matrix not found in folder : " << outFolder << std::endl;
	return 1;
      }
    }
    lChi2Fit.initialiseLeastSquareFit();
  }
  
  //loop on events
  HGCSSEvent * event = 0;
  std::vector<HGCSSSamplingSection> * ssvec = 0;
  std::vector<HGCSSSimHit> * simhitvec = 0;
  std::vector<HGCSSRecoHit> * rechitvec = 0;
  std::vector<HGCSSGenParticle> * genvec = 0;
  unsigned nPuVtx = 0;

  lSimTree->SetBranchAddress("HGCSSEvent",&event);
  lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
  lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
  lSimTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);

  lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  if (lRecTree->GetBranch("nPuVtx")) lRecTree->SetBranchAddress("nPuVtx",&nPuVtx);
  unsigned nTooFar = 0;
  unsigned nNoCluster = 0;

  const unsigned nRemove = 12;
  std::vector<unsigned> lToRemove;
  //unsigned list[nRemove] = {25,27,15,1,10,3,18,5,12,7,23,20};
  //unsigned list[nRemove] = {29,28,27,26,25,24,23,22,21,20,19,18};
  unsigned list[nRemove] = {0,1,2,3,4,5,6,7,8,9,10,11};

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;

    lSimTree->GetEntry(ievt);
    lRecTree->GetEntry(ievt);
    if (dofit) {
      bool found = lChi2Fit.setTruthInfo(genvec,1);
      if (!found) continue;
      //mask layers in turn
      lToRemove.clear();
      for (unsigned r(0); r<nRemove+1;++r){
	bool good = true;
	if (!doMatrix) {
	  good = lChi2Fit.getInitialPosition(ievt,nPuVtx,rechitvec,nTooFar,nNoCluster,lToRemove);
	  if (good && r==0) lChi2Fit.getOutTree()->Fill();
	}
	if (!good) continue;
	FitResult fit;
	//std::cout << r << " ";
	//fit.Print();
	if ( lChi2Fit.performLeastSquareFit(ievt,fit,lToRemove)==0 ){
	  //SignalEnergy.fillEnergies(ievt,(*ssvec),(*simhitvec),(*rechitvec),nPuVtx,fit);
	}
	else std::cout << " -- remove " << r << " Fit failed." << std::endl;
	if (r<nRemove) lToRemove.push_back(list[r]);
      }

    }
    else SignalEnergy.fillEnergies(ievt,(*ssvec),(*simhitvec),(*rechitvec),nPuVtx);

  }//loop on entries

  //finalise

  if (dofit) lChi2Fit.finaliseFit();
  SignalEnergy.finalise();

  outputFile->Write();
  //outputFile->Close();
  if (!doMatrix) std::cout << " -- Number of events with closest cluster 1 away from truth : " << nTooFar << std::endl;
  if (!doMatrix) std::cout << " -- Number of events with no cluster : " << nNoCluster << std::endl;
  std::cout << " - End of egammaResolution program." << std::endl;

  return 0;
  

}//main
