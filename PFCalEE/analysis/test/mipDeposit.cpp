#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TStyle.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"

#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSInfo.hh"
#include "HGCSSGeometryConversion.hh"
#include "HGCSSDetector.hh"


#include "utilities.h"

int main(int argc, char** argv){//main

  const unsigned nPar = static_cast<unsigned>(argc);
  if (nPar < 9) {
    std::cout << " Usage: "
              << argv[0] << " <nEvts to process (0=all)>"
	      << " <eta value: 1.600,2.000,2.850>"
	      << " <et value: 0, 1>"
	      << " <input eos path>"
	      << " <output dir>"
	      << " <BON ? 1>" 
	      << " <useDigi>"
	      << " <nRuns>"
	      << std::endl;
    return 1;
  }
  
  unsigned nEvts = atoi(argv[1]);
  std::string etaval = argv[2];
  std::string etval = argv[3];
  std::string inFilePath = argv[4];
  std::string plotBase = argv[5];
  bool bfieldON = atoi(argv[6]);
  bool useDigi = atoi(argv[7]);
  unsigned nRuns = atoi(argv[8]);


  //std::string inFilePath = "root://eoscms//eos/cms/store/group/dpg_hgcal/comm_hgcal/amagnan/HGCalTDR/gitV08-06-00/mu-/";
  //if (!useDigi) inFilePath += "HGcal_";
  //else inFilePath += "DigiIC3_eta1.90_2.10_200u";
  inFilePath += "_version63_model2_";
  if (bfieldON) inFilePath += "BON";
  else inFilePath += "BOFF";
  //if (etaval=="1.600") inFilePath += "_et20_eta";
  //else inFilePath += "_et13_eta";
  inFilePath += "_et";

  double expMPV = etaval=="1.600" ? 0.09 : etaval=="2.000"? 0.055 : 0.03;

  //const unsigned nSi = etaval=="1.600" ? 3 : etaval=="2.000"? 2 : 1;
  const unsigned firstLayerScint = 52;

  std::string outFilePath = plotBase+"/mipcalib_Et"+etval+"GeV_eta"+etaval;
  if (bfieldON) outFilePath += "BON";
  if (useDigi) outFilePath += "_fromDigi";
  outFilePath += ".root";

  TChain *lTree = new TChain(useDigi?"RecoTree":"HGCSSTree");
  TFile * mipFile = 0;
  for (unsigned i(0);i<nRuns;++i){
    std::ostringstream lstrrec;
    lstrrec << inFilePath << etval << "_eta" << etaval << "_run" << i << ".root";
    if (!testInputFile(lstrrec.str(),mipFile)) continue;
    lTree->AddFile(lstrrec.str().c_str());
  }
  
  
  if (!lTree){
    std::cout << " -- Error, tree MipTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  
  std::cout << " Trees added." << std::endl;
  
  unsigned versionNumber = 63;
  unsigned model = 2;
  unsigned shape = 1;
  double cellSize = 6.496345;
  double calorSizeXY = 5600;
  if (!useDigi){
    HGCSSInfo * info=(HGCSSInfo*)mipFile->Get("Info");

    assert(info);
    versionNumber = info->version();
    model = info->model();
    shape = info->shape();
    cellSize = info->cellSize();
    calorSizeXY = info->calorSizeXY();
    
  }

  HGCSSDetector & myDetector = theDetector();
  std::cout << " -- Calor size XY = " << calorSizeXY
	    << ", version number = " << versionNumber 
	    << ", model = " << model
	    << ", cellSize = " << cellSize
	    << ", shape = " << shape
	    << std::endl;

  bool bypassR = false;
  myDetector.buildDetector(versionNumber,model,true,false,bypassR);

  //initialise calibration class
  HGCSSCalibration mycalib(inFilePath,bypassR,2);
  HGCSSGeometryConversion geomConv(model,cellSize,bypassR,2);
  geomConv.setXYwidth(calorSizeXY);
  geomConv.setVersion(versionNumber);

  if (shape==2) geomConv.initialiseDiamondMap(calorSizeXY,10.);
  else if (shape==3) geomConv.initialiseTriangleMap(calorSizeXY,10.*sqrt(2.));
  else if (shape==1) geomConv.initialiseHoneyComb(calorSizeXY,cellSize);
  else if (shape==4) geomConv.initialiseSquareMap(calorSizeXY,10.);
  geomConv.initialiseSquareMap1(1.3,3.0,-1.*TMath::Pi(),TMath::Pi(),TMath::Pi()*2./360.);//eta phi segmentation, hardcoded '1.3' to include outer edge fo BH
  geomConv.initialiseSquareMap2(1.3,3.0,-1.*TMath::Pi(),TMath::Pi(),TMath::Pi()*2./288.);//eta phi segmentation, hardcoded '1.3' to include outer edge fo BH

  
  TFile *outputFile = TFile::Open(outFilePath.c_str(),"RECREATE");
  
  if (!outputFile) {
    std::cout << " -- Error, output file " << outFilePath << " cannot be opened. Please create output directory. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
  }
  outputFile->cd();
  
  const unsigned nLayers = 69;//33;
  
  TTree *outtree = new TTree("MipFit","Tree for Mip calibration");
  std::vector<double> EperLayer;
  EperLayer.resize(nLayers,0);
  std::vector<double> nHitsperLayer;
  nHitsperLayer.resize(nLayers,0);
  for (unsigned iL(0); iL<nLayers; ++iL){
    EperLayer[iL] = 0;
    std::ostringstream label;
    label << "E_layer" << iL;
    outtree->Branch(label.str().c_str(),&EperLayer[iL]);
    nHitsperLayer[iL] = 0;
    label.str("");
    label << "nHits_layer" << iL;
    outtree->Branch(label.str().c_str(),&nHitsperLayer[iL]);    
  }

  TH2F *p_nHits = new TH2F("nHits","; layer; Number of hits; Events",nLayers,0,nLayers,20,0,20);
  TH1F *p_hitEnergy_si = new TH1F("hitEnergy_si",";E (MeV);SimHits",1000,0,0.2);
  TH1F *p_hitEnergySel_si = new TH1F("hitEnergySel_si",";E (MeV);SimHits",1000,0,0.2);
  TH1F *p_hitEnergySel_si_EE = new TH1F("hitEnergySel_si_EE",";E (MeV);SimHits",1000,0,0.2);
  TH1F *p_hitEnergy_scint = new TH1F("hitEnergy_scint",";E (MeV);SimHits",1000,0,10);
  TH1F *p_hitEnergySel_scint = new TH1F("hitEnergySel_scint",";E (MeV);SimHits",1000,0,10);

  std::vector<HGCSSSimHit> * simhitvec = 0;
  if (!useDigi) lTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
  std::vector<HGCSSRecoHit> * rechitvec = 0;
  if (useDigi) lTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);

  if (nEvts==0 || nEvts > lTree->GetEntries()) nEvts = lTree->GetEntries();

  std::cout << " -- Processing " << nEvts << " entries." << std::endl;
  
  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    
    lTree->GetEntry(ievt);
    
    if (ievt%1000==0) std::cout << "- Event " << ievt << std::endl;
    
    for (unsigned iL(0); iL<nLayers; ++iL){
      nHitsperLayer[iL] = 0;
      EperLayer[iL] = 0;
    }
    
    unsigned nH = useDigi ? (*rechitvec).size() : (*simhitvec).size();

    std::vector<bool> passSel;
    passSel.resize(nH,false);

    for (unsigned iH(0); iH<nH; ++iH){//loop on hits
      double energy = 0;
      unsigned layer = 0;
      if (!useDigi) {
	HGCSSSimHit lHit = (*simhitvec)[iH];
	energy = lHit.energy();
	if (energy<=0) continue;
	layer = lHit.layer();
	const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
	DetectorEnum type = subdet.type;
	std::pair<double,double> xy = lHit.get_xy(subdet,geomConv,shape);
	double posx = xy.first;//lHit.get_x(cellSize);
	double posy = xy.second;//lHit.get_y(cellSize);
	double radius = sqrt(pow(posx,2)+pow(posy,2));
	unsigned nSi = 3;//geomConv.getNumberOfSiLayers(type,radius);//nSi
	if (lHit.silayer() >= nSi) continue;
      } else {
	HGCSSRecoHit lHit = (*rechitvec)[iH];
	energy = lHit.energy();
	if (energy<=0) continue;
	layer = lHit.layer();
      }
      if (layer < firstLayerScint) p_hitEnergy_si->Fill(energy);
      else p_hitEnergy_scint->Fill(energy);
      nHitsperLayer[layer]++;
      passSel[iH] = true;
    }//loop on hits

    for (unsigned iL(0); iL<nLayers; ++iL){
      p_nHits->Fill(iL,nHitsperLayer[iL]);
    }
    
    for (unsigned iH(0); iH<nH; ++iH){//loop on hits
      if (!passSel[iH]) continue;
      double energy = 0;
      unsigned layer = 0;
      if (!useDigi){
	HGCSSSimHit lHit = (*simhitvec)[iH];
	energy = lHit.energy();
	layer = lHit.layer();
	//if (layer < firstLayerScint && nHitsperLayer[layer]!=geomConv.getNumberOfSiLayers(type,radius)) continue;
	//if (layer > firstLayerScint && nHitsperLayer[layer]!=1) continue;
      } else {
	HGCSSRecoHit lHit = (*rechitvec)[iH];
	energy = lHit.energy();
	layer = lHit.layer();
      }
      EperLayer[layer] += energy;
      //std::cout << " -- hit " << iH << " layer " << lHit.layer() << " siLayer " << lHit.silayer() << " E " << lHit.energy() << " cellId " << lHit.cellid() << " Elayer = " <<  EperLayer[layer] << std::endl;
    }

    for (unsigned iL(0); iL<nLayers; ++iL){
      if (EperLayer[iL] <0.001) continue;
      if (iL < firstLayerScint) p_hitEnergySel_si->Fill(EperLayer[iL]);
      if (iL < 28) p_hitEnergySel_si_EE->Fill(EperLayer[iL]);
      else  p_hitEnergySel_scint->Fill(EperLayer[iL]);
    }

    outtree->Fill();
  }//loop on entries
  

  TCanvas *myc = new TCanvas("myc","myc",1);
  
  myc->cd();
  gPad->SetLogz(1);
  gStyle->SetOptStat(1111110);
  p_nHits->Draw("colz");

  myc->Update();
  myc->Print((plotBase+"/mipHits.pdf").c_str());
  myc->Print((plotBase+"/mipHits.C").c_str());


  myc->cd();
  gPad->SetLogy(1);
  gStyle->SetOptStat(1111110);
  gStyle->SetOptFit(1111);
  p_hitEnergy_si->Draw();
  p_hitEnergy_si->Fit("landau","+","",expMPV*0.8,expMPV*2);

  myc->Update();
  myc->Print((plotBase+"/mipDepositAll_si.pdf").c_str());
  myc->Print((plotBase+"/mipDepositAll_si.C").c_str());

  myc->cd();
  gPad->SetLogy(1);
  gStyle->SetOptStat(1111110);
  gStyle->SetOptFit(1111);
  p_hitEnergySel_si->Draw();
  p_hitEnergySel_si->Fit("landau","L+","",expMPV*0.8,expMPV*2);

  myc->Update();
  myc->Print((plotBase+"/mipDepositSel_si.pdf").c_str());
  myc->Print((plotBase+"/mipDepositSel_si.C").c_str());

  myc->cd();
  gPad->SetLogy(1);
  gStyle->SetOptStat(1111110);
  gStyle->SetOptFit(1111);
  p_hitEnergySel_si_EE->Draw();
  p_hitEnergySel_si_EE->Fit("landau","L+","",expMPV*0.8,expMPV*2);

  myc->Update();
  myc->Print((plotBase+"/mipDepositSel_si_EE.pdf").c_str());
  myc->Print((plotBase+"/mipDepositSel_si_EE.C").c_str());


  myc->cd();
  gPad->SetLogy(1);
  gStyle->SetOptStat(1111110);
  gStyle->SetOptFit(1111);
  p_hitEnergy_scint->Draw();
  p_hitEnergy_scint->Fit("landau","+");
  
  myc->Update();
  myc->Print((plotBase+"/mipDepositAll_scint.pdf").c_str());
  myc->Print((plotBase+"/mipDepositAll_scint.C").c_str());

  myc->cd();
  gPad->SetLogy(1);
  gStyle->SetOptStat(1111110);
  gStyle->SetOptFit(1111);
  p_hitEnergySel_scint->Draw();
  p_hitEnergySel_scint->Fit("landau","L+");
  
  myc->Update();
  myc->Print((plotBase+"/mipDepositSel_scint.pdf").c_str());
  myc->Print((plotBase+"/mipDepositSel_scint.C").c_str());

  outputFile->cd();
  outtree->Write();
  outputFile->Write();
  outputFile->Close();

  return 0;

}//main
