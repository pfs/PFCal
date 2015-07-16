#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#include <boost/algorithm/string.hpp>
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "boost/function.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"
#include "TGaxis.h"
//#include "macros/TDRStyle.h"
#include "TPaveStats.h"
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

#include "utilities.h"

using boost::lexical_cast;
namespace po=boost::program_options;

int main(int argc, char** argv){//main

  if (argc<8){
    std::cout << " Usage: "
	      << argv[0] 
	      << " <mipCut-200um>"<< std::endl
	      << " <eosPath>"<< std::endl
	      << " <version>"<< std::endl
	      << " <pt>"<< std::endl
	      << " <eta>"<< std::endl
	      << " <nRuns>" << std::endl
	      << " <nEvts to process (0=all)>"<< std::endl
	      << std::endl;
    return 1;
  }

  double mipCut = atoi(argv[1]);//200um Si: 26
  std::string eosPath = argv[2];
  unsigned version = 0;
  double pt = 0;
  std::istringstream(argv[3])>>version;
  std::istringstream(argv[4])>>pt;
  std::string eta = argv[5];
  unsigned nRuns = 0;
  std::istringstream(argv[6])>>nRuns;
  const unsigned pNevts = atoi(argv[7]);

  std::ostringstream inFilePath;
  inFilePath << eosPath << "HGcal_fixedVtx_version" << version << "_model2_BOFF_et" << pt << "_eta" << eta;
  if (eta=="1.75" || eta=="2.25" || eta=="2.75") inFilePath << "0";
  else inFilePath << "00";
  if (nRuns==0) inFilePath << ".root";
  else inFilePath << "_run0.root";

  std::ostringstream inputrec;
  inputrec << eosPath << "DigiIC3_fixedVtx_version" << version << "_model2_BOFF_et" << pt << "_eta" << eta;
  if (eta=="1.75" || eta=="2.25" || eta=="2.75") inputrec << "0";
  else inputrec << "00";

  std::ostringstream output;
  output << "NAbove" << mipCut << ".root";

  TFile *fout = TFile::Open(output.str().c_str(),"RECREATE");
  fout->cd();
  TTree *outtree = new TTree("outtree","Output tree with cells above 26 mips");
  unsigned nover = 0;
  unsigned ntot = 0;
  outtree->Branch("nover",&nover);
  outtree->Branch("ntot",&ntot);
  TH1F *nAbove = new TH1F("nAbove",";n_{cells}(E>26 MIPs);Probability",20,0,20);
  TH1F *nTot = new TH1F("nTot",";n_{cells};events",200,0,1000);
  TH2D *hxy = new TH2D("hxy",";x (mm);y (mm);cells E>26 MIPs",339,-1695,1695,339,-1695,1695);



  //SetTdrStyle();
  gStyle->SetPadRightMargin(0.1);

  //TFile *input = TFile::Open("/afs/cern.ch/work/a/amagnan/public/HGCalGeant4/git_V00-03-03/version_33/model_2/pi+/BON/et_3/a_2.000/run_0/DigiPFcal.root");
  //TFile *input = TFile::Open("root://eoscms//eos/cms/store/cmst3/group/hgcal/HGCalGeant4/gitV00-03-07/pi+/DigiIC2_version33_model2_BON_et2_eta2.000.root");
  //TFile *input = TFile::Open("root://eoscms//eos/cms/store/cmst3/group/hgcal/HGCalEEGeant4/gitV00-02-14/Hgg/DigiIC2_version12_model2_BOFF_run0.root");
  TChain *lRecTree = new TChain("RecoTree");


  TFile *inputSim = TFile::Open(inFilePath.str().c_str());
  if (!inputSim) {
    std::cout << " -- Error, input file " << inFilePath.str() << " not found. Exiting..." << std::endl;
    return 1;
  }
  HGCSSInfo* info =(HGCSSInfo*)inputSim->Get("Info");
  const unsigned versionNumber = info->version();
  //initialise detector
  HGCSSDetector & myDetector = theDetector();
 
  myDetector.buildDetector(versionNumber,true,false);
  HGCSSGeometryConversion geomConv(inFilePath.str(),info->model(),info->cellSize());

  const unsigned nLayers = myDetector.nLayers();

  if (nRuns==0) {
    lRecTree->AddFile((inputrec.str()+".root").c_str());
  } else {
    for (unsigned i(0);i<nRuns;++i){
      std::ostringstream lstrrec;
      lstrrec << inputrec.str() << "_run" << i << ".root";
      lRecTree->AddFile(lstrrec.str().c_str());
    }
  }

  gStyle->SetOptStat("eMRuo");

  const unsigned nEvts = (pNevts > lRecTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lRecTree->GetEntries()) : pNevts;

  std::cout << " -- Processing " << nEvts << " events." << std::endl
	    << " -- mipCut is set to " << mipCut << " mips for 200um Si." << std::endl;

  TCanvas *myc = new TCanvas("myc","myc",1500,1000);
  TCanvas *myc2 = new TCanvas("myc2","myc2",1500,1000);
  std::ostringstream lvar;
  lvar << "HGCSSRecoHitVec.ypos_:HGCSSRecoHitVec.xpos_>>hxy";
  std::ostringstream lcut;
  lcut << "HGCSSRecoHitVec.energy_>" << mipCut << " && HGCSSRecoHitVec.layer_<" << nLayers-12;

  myc->cd();
  hxy->GetZaxis()->SetTitleOffset(0.5);
  hxy->SetStats(0);
  lRecTree->Draw(lvar.str().c_str(),lcut.str().c_str());

  hxy->Draw("colz");
  myc->Update();
  std::ostringstream lsave;
  lsave << "yvsx_pi_Eabove" << mipCut << ".pdf";
  myc->Print(lsave.str().c_str());



  std::vector<HGCSSRecoHit> * rechitvec = 0;
  lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (ievt%100 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
    
    lRecTree->GetEntry(ievt);

    nover = 0;
    ntot = 0;

    for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
      const HGCSSRecoHit & lHit = (*rechitvec)[iH];
      double energy = lHit.energy();
      double leta = lHit.eta();
      //if (leta>2.1 || leta < 1.9) continue;
      unsigned layer = lHit.layer();
      const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
      if (subdet.type == DetectorEnum::BHCAL1 || subdet.type == DetectorEnum::BHCAL2) continue;
      double lRadius = sqrt(pow(lHit.get_x(),2)+pow(lHit.get_y(),2));
      double cutVal = mipCut*2./geomConv.getNumberOfSiLayers(subdet.type,lRadius);
      if (energy>cutVal){
	nover++;
      }
      ntot++;
    }//loop on hits
    nAbove->Fill(nover);
    nTot->Fill(ntot);

    outtree->Fill();

  }//loop on entries
  
  nAbove->Scale(1./nEvts);
  myc2->cd();
  //gPad->SetGridx(1);
  //gPad->SetGridy(1);
  nAbove->GetXaxis()->SetRangeUser(0,200);
  gPad->SetLogy(1);
  nAbove->Draw("histtext");

  TLatex lat;
  lsave.str("");
  lsave << "#pi^{-} p_{T}=" << pt << " GeV, #eta = " << eta;
  lat.DrawLatexNDC(0.3,0.95,lsave.str().c_str());
  //lat.DrawLatexNDC(0.3,0.95,"Pythia H#rightarrow#gamma#gamma 1.9 < #eta < 2.1");
  lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 Standalone Simulation");
  myc2->Update();
  lsave.str("");
  lsave << "nAbove" << mipCut << "mips_pi" << pt << ".pdf";

  myc2->Print(lsave.str().c_str());

  fout->Write();
  return 0;

}//main
