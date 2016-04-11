#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#include<cmath>

#include "TFile.h"
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

#include "TDRStyle.h"

double etafunc(double x, double y, double z) {
  double theta = acos(fabs(z)/sqrt(z*z+x*x+y*y));
  double leta = -log(tan(theta/2.));
  if (z>0) return leta;
  else return -leta;
}

int plotMomentumMuonsBHCAL(){


  const unsigned nFiles = 10;
  const unsigned nE = 4;
  const unsigned et[nE] = {5,10,20,50};
  const double eta=2.8;
  TCanvas *myc = new TCanvas("myc","myc",1);

  TLatex lat;
  char buf[500];

  for (unsigned iE(0); iE<nE; ++iE){//loop on energy points
    std::ostringstream lstr;
    lstr << "PLOTS/output_layer51_eta" << eta << "_et" << et[iE];
    lstr << ".root";
    TFile *outputFile = TFile::Open(lstr.str().c_str(),"RECREATE");
    
    if (!outputFile) {
      std::cout << " -- Error, output file " << lstr.str() << " cannot be opened. Please create output directory. Exiting..." << std::endl;
      return 1;
    }
    else {
      std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
    }
    outputFile->cd();

    TH1F *p_logmomentum = new TH1F("p_logmomentum",";log(p) (log(GeV));particles",100,1,3);
    TH1F *p_momentum  = new TH1F("p_momentum",";p (GeV);particles",5000,0,500);
    
    TH1F *p_trkphi = new TH1F("p_trkphi",";#phi;muons",1000,1,2);
    TH1F *p_phi = new TH1F("p_phi",";#phi;muons",1000,1,2);
    TH1F *p_eta = new TH1F("p_eta",";#eta;muons",100,2,3);
    TH1F *p_dphi = new TH1F("p_dphi",";#Delta#phi;muons",200,-0.2,0.4);
    TH1F *p_dphi_pos = new TH1F("p_dphi_pos",";atan(y/x);muons",200,1.5,1.7);
    TH1F *p_dphi_mom = new TH1F("p_dphi_mom",";atan(py/px);muons",200,1.5,1.9);
    TH2F *p_dphi_pos_vs_mom = new TH2F("p_dphi_pos_vs_mom",";atan(py/px);atan(y/x);muons",200,1.5,1.9,200,1.5,1.7);

    unsigned nEvts = 0;
    
    for (unsigned iF(0); iF<nFiles;++iF){
      
      std::ifstream input;
      std::ostringstream label;
      //label << "/afs/cern.ch/work/a/amagnan/public/HGCalGeant4/git_V00-03-03/version_25/model_2/minbias/BON/run_0/momentum_list_layer30.dat";
      label << "/afs/cern.ch/work/a/amagnan/public/HGCalMuons/git_V06b-03-05/version_33/model_2/mu-/BON//et_" << et[iE] << "/eta_" << eta << "00/run_" << iF << "/momentum_list_layer51.dat";
      input.open(label.str().c_str());
      if (!input.is_open()){
	std::cout << " Cannot open input file " << label.str() << ". Continue..." << std::endl;
	continue;
      }
      else std::cout << " Successfully opened file " << label.str() << std::endl;
      
      unsigned evtid=1000;
      int trkid=0;
      int pdgid=0;
      double time=0;
      double x=0;
      double y=0;
      double z=0;
      double px=0;
      double py=0;
      double pz=0;
      std::string volStart;
      std::string volEnd;
      
      while (!input.eof()){
	std::string event("");
	std::getline(input,event);
	//std::cout << " -- line: " << event << std::endl;
	if (event.find("Event")!=event.npos){
	  std::istringstream(event.substr(6,event.npos-6))>>evtid;
	  //std::cout << " --Processing event id " << evtid 
	  //<< " " << event.substr(6,event.npos-6)
	  //<< std::endl;
	  //return 1;
	  nEvts++;
	  continue;
	}
	std::istringstream iss(event);
	iss>>volStart>>volEnd>>trkid>>pdgid>>time>>x>>y>>z>>px>>py>>pz;
	double p = sqrt(px*px+py*py+pz*pz);
	p_momentum->Fill(p/1000.);
	p_logmomentum->Fill(log10(p/1000.));
	p_phi->Fill(std::atan2(y,x));
	p_trkphi->Fill(std::atan2(py,px));
	p_eta->Fill(etafunc(x,y,z));
	p_dphi->Fill(std::atan2(py,px)-std::atan2(y,x));
	p_dphi_pos->Fill(std::atan2(y,x));
	p_dphi_mom->Fill(std::atan2(py,px));
 	p_dphi_pos_vs_mom->Fill(std::atan2(py,px),std::atan2(y,x));
     }//read input file
      
    }//loop on files
    std::cout << " -- Read " << nEvts << " events." << std::endl;

    sprintf(buf,"#mu p_{T} = %d GeV, #eta=%3.1f",et[iE],eta);

    myc->cd();
    gStyle->SetOptStat("eMR");
    p_phi->Draw();
    lat.DrawLatexNDC(0.2,0.85,buf);
    std::ostringstream lsave;
    lsave << "PLOTS/phi_eta" << eta << "_et" << et[iE] << ".pdf";
    myc->Print(lsave.str().c_str());

    myc->cd();
    gStyle->SetOptStat("eMR");
    p_trkphi->Draw();
    lat.DrawLatexNDC(0.2,0.85,buf);
    lsave.str("");
    lsave << "PLOTS/trkphi_eta" << eta << "_et" << et[iE] << ".pdf";
    myc->Print(lsave.str().c_str());

    myc->cd();
    gStyle->SetOptStat("eMR");
    p_dphi->Draw();
    lat.DrawLatexNDC(0.2,0.85,buf);
    lsave.str("");
    lsave << "PLOTS/deltaphi_eta" << eta << "_et" << et[iE] << ".pdf";
    myc->Print(lsave.str().c_str());

    myc->cd();
    gStyle->SetOptStat("eMR");
    p_dphi_pos->Draw();
    lat.DrawLatexNDC(0.2,0.85,buf);
    lsave.str("");
    lsave << "PLOTS/deltaphi_position_eta" << eta << "_et" << et[iE] << ".pdf";
    myc->Print(lsave.str().c_str());

    myc->cd();
    gStyle->SetOptStat("eMR");
    p_dphi_mom->Draw();
    lat.DrawLatexNDC(0.2,0.85,buf);
    lsave.str("");
    lsave << "PLOTS/deltaphi_momentum_eta" << eta << "_et" << et[iE] << ".pdf";
    myc->Print(lsave.str().c_str());

    myc->cd();
    gStyle->SetOptStat("eMR");
    p_dphi_pos_vs_mom->Draw("colz");
    lat.DrawLatexNDC(0.2,0.85,buf);
    lsave.str("");
    lsave << "PLOTS/deltaphi_posvsmom_eta" << eta << "_et" << et[iE] << ".pdf";
    myc->Print(lsave.str().c_str());

    myc->cd();
    gStyle->SetOptStat("eMR");
    p_eta->Draw();
    lat.DrawLatexNDC(0.2,0.85,buf);
    lsave.str("");
    lsave << "PLOTS/eta_eta" << eta << "_et" << et[iE] << ".pdf";
    myc->Print(lsave.str().c_str());
  
    myc->cd();
    gStyle->SetOptStat("eMR");
    p_momentum->Draw();
    lat.DrawLatexNDC(0.2,0.85,buf);
    lsave.str("");
    lsave << "PLOTS/p_eta" << eta << "_et" << et[iE] << ".pdf";
    myc->Print(lsave.str().c_str());
  
    outputFile->Write();
  }//loop on energies
  return 0;
  
}//main
