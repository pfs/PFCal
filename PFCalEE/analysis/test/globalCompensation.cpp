#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TRandom3.h"

#include "HGCSSInfo.hh"
#include "HGCSSEvent.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSParameters.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDigitisation.hh"

#include "TRandom3.h"

double getResolution(unsigned genEn,unsigned versionNumber){
  double s = 0.73;
  double c = 0.044;
  double n = 0.02;
  if (versionNumber==23){
    s = 0.48;
    c = 0.045;
    n = 0.18;
  }
  return genEn*sqrt(pow(s/sqrt(genEn),2)+pow(n/genEn,2)+pow(c,2));
};

double showerEnergy(const double Cglobal,
		    const double Etotcal,
		    const bool correctLinearity){
  double a1 = 0.89;
  double a2 = 0.001;
  double a3 = 0;//1.2e-5;
  double Esh = Cglobal*Etotcal;
  if (!correctLinearity) return Esh;
  return Esh*(a1+a2*Esh+a3*Esh*Esh);
};

int main(int argc, char** argv){//main  

  if (argc < 7) {
    std::cout << " Usage: " 
	      << argv[0] << " <nEvts to process (0=all)>"
	      << " <path to input files>"
	      << " <name of input sim file>"
	      << " <name of input reco file>"
	      << " <full path to output file>"
	      << " <number of si layers to consider: 1,2 or 3>" 
	      << " <optional: debug (default=0)>"
	      << std::endl;
    return 1;
  }

  //////////////////////////////////////////////////////////
  //// Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////

  bool concept = true;

  bool selectEarlyDecay = true;

  const unsigned nLimits = 10;//5;
  const double pElim[nLimits] = {5,6,7,8,9,10,11,12,13,14};
  const unsigned idxRef = 0;

  double FHcalEMCalib = 40.4;//39.81;//38;
  double FHcalEMOffset = -3.9;//1.9;//-15;
  double BHcalEMCalib = FHcalEMCalib;//39.81;//38;
  double BHcalEMOffset = 0;//1.9;//-15;
  double HcalPionCalib = 1/1.19;
  double HcalPionOffset = -0.81;
  double BHcalSlope = 1;

  //////////////////////////////////////////////////////////
  //// End Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////

  TRandom3 *lRndm = new TRandom3();
  lRndm->SetSeed(1234);

  std::cout << " -- Random3 seed = " << lRndm->GetSeed() << std::endl
	    << " ----------------------------------------" << std::endl;

  //for basic control plots
  TCanvas *mycE = new TCanvas("mycE","mycE",1500,1000);
  gStyle->SetOptStat("eMRuoi");

  const unsigned pNevts = atoi(argv[1]);
  std::string filePath = argv[2];
  std::string simFileName = argv[3];
  std::string recoFileName = argv[4];

  std::string inFilePath = filePath+simFileName;

  std::string outPath = argv[5];
  unsigned nSiLayers = 2;
  nSiLayers = atoi(argv[6]);

  unsigned debug = 0;
  if (argc >7) debug = atoi(argv[7]);


  bool isEM = false;

  if (inFilePath.find("e-")!=inFilePath.npos || 
      inFilePath.find("e+")!=inFilePath.npos) isEM = true;


  std::size_t begin = inFilePath.find_last_of("_e")+1;
  std::size_t end = inFilePath.find(".root");
  unsigned genEn = 0;
  std::cout << inFilePath << " " << begin << " " << end << " " << inFilePath.substr(begin,end-begin) << std::endl;
  std::istringstream(inFilePath.substr(begin,end-begin))>>genEn;

  if (selectEarlyDecay && isEM) {
    selectEarlyDecay = false;
    HcalPionCalib = 1;
    HcalPionOffset = 0;
  }

  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << filePath << std::endl
	    << " -- Output file path: " << outPath << std::endl
	    << " -- Input energy is: " << genEn << " GeV." << std::endl
	    << " -- Requiring " << nSiLayers << " si layers." << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events." << std::endl;




  /////////////////////////////////////////////////////////////
  //input
  /////////////////////////////////////////////////////////////


  std::ostringstream input;
  input << filePath << "/" << simFileName;

  TFile *simFile = TFile::Open(input.str().c_str());

  if (!simFile) {
    std::cout << " -- Error, input file " << input.str() << " cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  else std::cout << " -- input file " << simFile->GetName() << " successfully opened." << std::endl;
  
  TTree *lSimTree = (TTree*)simFile->Get("HGCSSTree");
  if (!lSimTree){
    std::cout << " -- Error, tree HGCSSTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  input.str("");
  input << filePath << "/" << recoFileName;
  
  TFile *recFile = TFile::Open(input.str().c_str());

  if (!recFile) {
    std::cout << " -- Error, input file " << input.str() << " cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  else std::cout << " -- input file " << recFile->GetName() << " successfully opened." << std::endl;

  TTree *lRecTree = (TTree*)recFile->Get("RecoTree");
  if (!lRecTree){
    std::cout << " -- Error, tree RecoTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }


  /////////////////////////////////////////////////////////////
  //Info
  /////////////////////////////////////////////////////////////

  HGCSSInfo * info=(HGCSSInfo*)simFile->Get("Info");
  const double cellSize = info->cellSize();
  const unsigned versionNumber = info->version();
  const unsigned model = info->model();
  
  //models 0,1 or 3.
  bool isTBsetup = (model != 2);
  bool isCaliceHcal = versionNumber==23;//inFilePath.find("version23")!=inFilePath.npos || inFilePath.find("version_23")!=inFilePath.npos;


  std::cout << " -- Version number is : " << versionNumber 
	    << ", model = " << model
	    << ", cellSize = " << cellSize
	    << std::endl;


  //initialise detector
  HGCSSDetector & myDetector = theDetector();

  unsigned indices[7] = {0,0,0,0,0,0,0};
  //fill layer indices
  if (versionNumber==22) {
    indices[4] = 0;
    indices[5] = 10;
    indices[6] = 10;
  }
  else if (isCaliceHcal) {
    indices[3] = 0;
    indices[4] = 38;
    indices[5] = 47;
    indices[6] = 54;
  }
  else if (versionNumber==21) {
    indices[3] = 0;
    indices[4] = 24;
    indices[5] = 34;
    indices[6] = 34;
  }
  else if (versionNumber < 20){
    indices[0] = 0;
    indices[1] = versionNumber==8?11:10;
    indices[2] = versionNumber==8?21:20;
    indices[3] = versionNumber==8?31:30;
    indices[4] = indices[3];
    indices[5] = indices[3];
    indices[6] = indices[3];
  }
  else {
    indices[0] = 0;
    indices[1] = 11;
    indices[2] = 21;
    indices[3] = 31;
    indices[4] = 55;
    indices[5] = 65;
    indices[6] = 65;
  }

  if (versionNumber == 23) {
    FHcalEMCalib = 40.4;
    FHcalEMOffset = -3.9;
    BHcalEMCalib = FHcalEMCalib;
    BHcalEMOffset = 0;//1.9;//-15;
    HcalPionCalib = 0.84;
    HcalPionOffset = 0;//-1.4;//-0.81;
    BHcalSlope = 1;
  }

  myDetector.buildDetector(indices,concept,isCaliceHcal);

  //initialise calibration class
  HGCSSCalibration mycalib(inFilePath);
  HGCSSDigitisation myDigitiser;
  myDigitiser.setRandomSeed(lRndm->GetSeed());

  const unsigned nLayers = myDetector.nLayers();
  const unsigned nSections = myDetector.nSections();

  std::cout << " -- N layers = " << nLayers << std::endl
	    << " -- N sections = " << nSections << std::endl;


  TFile *outputFile = TFile::Open(outPath.c_str(),"RECREATE");
  
  if (!outputFile) {
    std::cout << " -- Error, output file " << outPath << " cannot be opened. Please create output directory. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
  }
  
  
  outputFile->cd();

  TH1F *p_ErecoTotal = new TH1F("p_ErecoTotal",";Ereco (GeV)",2000,0,1000);
  p_ErecoTotal->StatOverflows();
  p_ErecoTotal->Sumw2();
  TH1F *p_EshowerCor = new TH1F("p_EshowerCor",";Ecor (GeV)",2000,0,1000);
  p_EshowerCor->StatOverflows();
  p_EshowerCor->Sumw2();
  
  TH1F *p_hitSpectrum_lowTail =  new TH1F("p_hitSpectrum_lowTail",";E^{hit}_{HCAL} (MIPs)",1000,0,500);
  p_hitSpectrum_lowTail->StatOverflows();
  p_hitSpectrum_lowTail->Sumw2();
  TH1F *p_hitSpectrum_highTail =  new TH1F("p_hitSpectrum_highTail",";E^{hit}_{HCAL} (MIPs)",1000,0,500);
  p_hitSpectrum_highTail->StatOverflows();
  p_hitSpectrum_highTail->Sumw2();
  TH1F *p_meanHitSpectrum_lowTail =  new TH1F("p_meanHitSpectrum_lowTail",";<E^{hit}_{HCAL}> (MIPs)",1000,0,200);
  p_meanHitSpectrum_lowTail->StatOverflows();
  p_meanHitSpectrum_lowTail->Sumw2();
  TH1F *p_meanHitSpectrum_highTail =  new TH1F("p_meanHitSpectrum_highTail",";<E^{hit}_{HCAL}> (MIPs)",1000,0,200);
  p_meanHitSpectrum_highTail->StatOverflows();
  p_meanHitSpectrum_highTail->Sumw2();
  TH1F *p_Cglobal[nLimits];
  TH1F *p_Eshower[nLimits];
  TH2F *p_EvsCglobal[nLimits];

  //TH2F * EmipHits = new TH2F("EmipHits",";x(mm);y(mm)",50,-500,500,50,-500,500);
  //TH2F * EmipHits = new TH2F("EmipHits",";x(mm);y(mm)",16,-240,240,16,-240,240);

  std::ostringstream lName;
  for (unsigned iLim(0); iLim<nLimits;++iLim){
    lName.str("");
    lName << "p_Cglobal_" << iLim;
    p_Cglobal[iLim] =  new TH1F(lName.str().c_str(),";C_{global}",1000,0,2);
      
    lName.str("");
    lName << "p_Eshower_" << iLim;
    p_Eshower[iLim] = new TH1F(lName.str().c_str(),";E_{shower} (GeV)",2000,0,1000);

    lName.str("");
    lName << "p_EvsCglobal_" << iLim;
    p_EvsCglobal[iLim] =  new TH2F(lName.str().c_str(),";C_{global};Etot (GeV)",
				   1000,0,2,
				   2000,0,1000);

  }

  HGCSSEvent * event = 0;
  std::vector<HGCSSSamplingSection> * ssvec = 0;
  std::vector<HGCSSSimHit> * simhitvec = 0;
  std::vector<HGCSSRecoHit> * rechitvec = 0;
  
  lSimTree->SetBranchAddress("HGCSSEvent",&event);
  lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
  lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
  
  lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);

  const unsigned nEvts = ((pNevts > lSimTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lSimTree->GetEntries()) : pNevts) ;
  
  std::cout << "- Processing = " << nEvts  << " events out of " << lSimTree->GetEntries() << std::endl;
  

  //double Esim[nSections];
  double Ereco[nSections];
  for (unsigned iD(0); iD<nSections; ++iD){
    //Esim[iD] = 0;
    Ereco[iD] = 0;
  }


  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
    
    lSimTree->GetEntry(ievt);
    lRecTree->GetEntry(ievt);

    if (debug){
      std::cout << "... Size of hit vectors: sim = " <<  (*simhitvec).size() << ", reco = " << (*rechitvec).size()<< std::endl;
    }

    unsigned firstInteraction = 0;
    for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
      HGCSSSimHit lHit = (*simhitvec)[iH];

      //discard some si layers...
      if (lHit.silayer() >= nSiLayers) continue; 

      unsigned layer = lHit.layer();
      if ( firstInteraction == 0 &&
	   (lHit.nNeutrons()>0 || 
	    lHit.nProtons()>0 ||
	    lHit.nHadrons()>0 ) && 
	   lHit.mainParentTrackID() > 0
	   ) firstInteraction = layer;
    }


    //get mean hit energy
    double EmipMean = 0;
    unsigned nHits = 0;
    
    for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on rechits
      HGCSSRecoHit lHit = (*rechitvec)[iH];
      if (debug>1) {
	std::cout << " --  RecoHit " << iH << "/" << (*rechitvec).size() << " --" << std::endl
		  << " --  position x,y " << lHit.get_x() << "," << lHit.get_y() << std::endl;
	lHit.Print(std::cout);
      }
      
      double posx = lHit.get_x();
      double posy = lHit.get_y();
      //double posz = lHit.get_z();
      
      double energy = lHit.energy();
      unsigned layer = lHit.layer();

      if (layer >= nLayers) {
	std::cout << " WARNING! RecoHits with layer " << layer << " outside of detector's definition range ! Please fix the digitiser or the detector definition used here. Ignoring..." << std::endl;
	continue;
      }
      unsigned sec =  myDetector.getSection(layer);
            
      //p_recoxy[layer]->Fill(posx,posy,energy);
      //EtotRec[layer] += energy;
      //if (debug>1) std::cout << "-hit" << iH << "-" << layer << " " << energy << " " << EtotRec[layer];

      double absweight = (*ssvec)[layer].volX0trans()/(*ssvec)[0].volX0trans();
      
      Ereco[sec] += energy*absweight;
      //if (sec<1) EmipHits->Fill(posx,posy,energy);

      if (sec==0){//do just for FHCAL
	EmipMean += energy;
	nHits++;
      }

    }//loop on rechits
    
    double Eecal = 0;
    if (myDetector.section(DetectorEnum::FECAL)<nSections) Eecal += myDigitiser.MIPtoGeV(myDetector.subDetectorByEnum(DetectorEnum::FECAL),Ereco[myDetector.section(DetectorEnum::FECAL)]);
    if (myDetector.section(DetectorEnum::MECAL)<nSections) Eecal += myDigitiser.MIPtoGeV(myDetector.subDetectorByEnum(DetectorEnum::MECAL),Ereco[myDetector.section(DetectorEnum::MECAL)]);
    if (myDetector.section(DetectorEnum::BECAL)<nSections) Eecal += myDigitiser.MIPtoGeV(myDetector.subDetectorByEnum(DetectorEnum::BECAL),Ereco[myDetector.section(DetectorEnum::BECAL)]);
    /*double Efhcal = 0;
    if (myDetector.section(DetectorEnum::FHCAL)<nSections) Efhcal += myDigitiser.MIPtoGeV(myDetector.subDetectorByEnum(DetectorEnum::FHCAL),Ereco[myDetector.section(DetectorEnum::FHCAL)]);
    double Ebhcal = 0;
    if (myDetector.section(DetectorEnum::BHCAL1)<nSections) Ebhcal += myDigitiser.MIPtoGeV(myDetector.subDetectorByEnum(DetectorEnum::BHCAL1),Ereco[myDetector.section(DetectorEnum::BHCAL1)]);
    if (myDetector.section(DetectorEnum::BHCAL2)<nSections) Ebhcal += myDigitiser.MIPtoGeV(myDetector.subDetectorByEnum(DetectorEnum::BHCAL2),Ereco[myDetector.section(DetectorEnum::BHCAL2)]);*/
    double Efhcal = 0;
    if (myDetector.section(DetectorEnum::FHCAL)<nSections) Efhcal += (Ereco[myDetector.section(DetectorEnum::FHCAL)]-FHcalEMOffset)/FHcalEMCalib;
    double Ebhcal = 0;
    if (myDetector.section(DetectorEnum::BHCAL1)<nSections) Ebhcal += (Ereco[myDetector.section(DetectorEnum::BHCAL1)]-BHcalEMOffset)/BHcalEMCalib;
    if (myDetector.section(DetectorEnum::BHCAL2)<nSections) Ebhcal += (Ereco[myDetector.section(DetectorEnum::BHCAL2)]-BHcalEMOffset)/BHcalEMCalib;

    double Etotcal = Eecal+(Efhcal+(1./BHcalSlope*Ebhcal)-HcalPionOffset)/HcalPionCalib;

    bool doFill = true;
    if (selectEarlyDecay && firstInteraction>5) doFill = false;

   if (doFill){
     //p_EsimTotal->Fill(etotmips);
      p_ErecoTotal->Fill(Etotcal);
    }
 
   if (doFill) {
     unsigned nHitsCountNum[nLimits];
     unsigned nHitsCountDen = 0;

     EmipMean = EmipMean/nHits;
     bool lowTail = Etotcal < (genEn-getResolution(genEn,versionNumber));
     bool highTail = Etotcal > (genEn+getResolution(genEn,versionNumber));
     //fill Cglobal histos
     if (lowTail) p_meanHitSpectrum_lowTail->Fill(EmipMean);
     if (highTail) p_meanHitSpectrum_highTail->Fill(EmipMean);

     for (unsigned iLim(0); iLim<nLimits;++iLim){
       nHitsCountNum[iLim] = 0;
       for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on rechits
	 HGCSSRecoHit lHit = (*rechitvec)[iH];
	 double energy = lHit.energy();
	 unsigned sec =  myDetector.getSection(lHit.layer());
	 if (sec==0){//use just fhcal
	   if (energy<=pElim[iLim]) nHitsCountNum[iLim]++;
	   //do just once
	   if (iLim==0){
	     if (energy<=EmipMean) nHitsCountDen++;
	     if (lowTail) {
	       p_hitSpectrum_lowTail->Fill(energy);
	     }
	     if (highTail) {
	       p_hitSpectrum_highTail->Fill(energy);
	     }
	   }
	 }
       }
       double Cglobal = 0;
       if (nHitsCountDen>0) Cglobal = nHitsCountNum[iLim]*1.0/nHitsCountDen;
       //if (iLim==idxRef && Cglobal>0.999 && Cglobal<1.001) std::cout << nHitsCountNum[iLim] << "/" << nHitsCountDen << " = " << Cglobal << std::endl;
       p_Cglobal[iLim]->Fill(Cglobal);
       p_EvsCglobal[iLim]->Fill(Cglobal,Efhcal);
       p_Eshower[iLim]->Fill(showerEnergy(Cglobal,Etotcal,false));
       if (iLim==idxRef) p_EshowerCor->Fill(showerEnergy(Cglobal,Etotcal,true));
     }
   }//if fill
   
   for (unsigned iD(0); iD<nSections; ++iD){
     Ereco[iD] = 0;
   }
   
   //EmipHits->Reset();
  }//loop on events


  std::cout << " -- Summary of total Energy: " << std::endl
	    << " ---- entries " << p_ErecoTotal->GetEntries() 
	    << " mean " << p_ErecoTotal->GetMean() 
	    << " rms " << p_ErecoTotal->GetRMS() 
	    << " underflows " << p_ErecoTotal->GetBinContent(0)
	    << " overflows " << p_ErecoTotal->GetBinContent(p_ErecoTotal->GetNbinsX()+1)
	      << std::endl
	      << " -- Summary of total Energy corrected: " << std::endl
	      << " ---- entries " << p_EshowerCor->GetEntries() 
	      << " mean " << p_EshowerCor->GetMean() 
	      << " rms " << p_EshowerCor->GetRMS() 
	      << " underflows " << p_EshowerCor->GetBinContent(0)
	      << " overflows " << p_EshowerCor->GetBinContent(p_EshowerCor->GetNbinsX()+1)
	      << std::endl
	      << " -- Summary of hit spectra low tail: " << std::endl
	      << " ---- entries " << p_hitSpectrum_lowTail->GetEntries() 
	      << " mean " << p_hitSpectrum_lowTail->GetMean() 
	      << " rms " << p_hitSpectrum_lowTail->GetRMS() 
	      << " underflows " << p_hitSpectrum_lowTail->GetBinContent(0)
	      << " overflows " << p_hitSpectrum_lowTail->GetBinContent(p_hitSpectrum_lowTail->GetNbinsX()+1)
	      << std::endl
	      << " -- Summary of hit spectra high tail: " << std::endl
	      << " ---- entries " << p_hitSpectrum_highTail->GetEntries() 
	      << " mean " << p_hitSpectrum_highTail->GetMean() 
	      << " rms " << p_hitSpectrum_highTail->GetRMS() 
	      << " underflows " << p_hitSpectrum_highTail->GetBinContent(0)
	      << " overflows " << p_hitSpectrum_highTail->GetBinContent(p_hitSpectrum_highTail->GetNbinsX()+1)
	      << std::endl;




  outputFile->Write();
  //outputFile->Close();
  
  return 0;
  
  
}//main
