#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>

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

int plotMomentumMinbias(){

  SetTdrStyle();


  TH1F *hP = 0;
  TH1F *hPi = 0;
  TH1F *hN = 0;
  TH1F *hSi = 0;

  unsigned siLayer=2;
  std::string silabel = "_si2";

  const unsigned nF = 2;
  TFile *f[nF];

  f[0] = TFile::Open("../PLOTS/output_hepmc.root");
  //f[1] = TFile::Open("../PLOTS/hadronMomentumAtShowerMax_all.root");

  std::ostringstream lstr;
  lstr << "PLOTS/output_layer30_";
  if (siLayer==3) lstr << "siAll";
  else lstr << "si" << siLayer;
  lstr << ".root";

  f[1] = TFile::Open(lstr.str().c_str());

  f[1]->cd();
  TGraph *gr = (TGraph*)gDirectory->Get("grID");
  f[1]->cd("neutrons");
  TH2F *hTvsLogp = (TH2F*)gDirectory->Get("p_timevslogp_n");
  TCanvas *myc[2*nF+2];

  myc[0] = new TCanvas("mycIP","IP",1);
  myc[1] = new TCanvas("mycSM","Shower Max",1);
  myc[2] = new TCanvas("mycIPlog","IP",1);
  myc[3] = new TCanvas("mycSMlog","Shower Max",1);
  myc[4] = new TCanvas("mycID","ID",1500,1000);
  myc[5] = new TCanvas("mycT","T",1000,500);

  gStyle->SetOptStat(0);
  TLatex lat;

  std::string label[2*nF] = {
    "hadronMomentumAtIP",
    "hadronMomentumAtShowerMax",
    "hadronMomentumAtIP_log",
    "hadronMomentumAtShowerMax_log"
  };

  for (unsigned iF(0);iF<nF;++iF){//loop on files
    if (!f[iF]) {
      std::cout << " File not found. Label: " << label[iF] << " Continue..." << std::endl;
      continue;
    }
    for (unsigned dolog(0); dolog<2;++dolog){
      f[iF]->cd();
      if (iF==0){
	hP = (TH1F*)gDirectory->Get(dolog==0?"hProton":"hProtonLog");
	hPi = (TH1F*)gDirectory->Get(dolog==0?"hPipm":"hPipmLog");
	hN = (TH1F*)gDirectory->Get(dolog==0?"hNeutron":"hNeutronLog");
	if (!hP || !hN || !hPi) {
	  std::cout << " -- histos not found ! Continue..." << std::endl;
	  continue;
	}
      }
      else {
	f[iF]->cd("protons");
	hP = (TH1F*)gDirectory->Get(dolog==0?"p_momentum_p":"p_logmomentum_p");
	f[iF]->cd();
	f[iF]->cd("neutrons");
	hN = (TH1F*)gDirectory->Get(dolog==0?"p_momentum_n":"p_logmomentum_n");
	f[iF]->cd();
	f[iF]->cd("pipm");
	hPi = (TH1F*)gDirectory->Get(dolog==0?"p_momentum_pi":"p_logmomentum_pi");
	f[iF]->cd();
	f[iF]->cd("Si");
	hSi = (TH1F*)gDirectory->Get(dolog==0?"p_momentum_Si":"p_logmomentum_Si");
	if (!hP || !hN || !hPi){
	  std::cout << " -- histos not found ! Continue..." << std::endl;
	  continue;
	}
      }
      myc[2*dolog+iF]->cd();
      myc[2*dolog+iF]->SetGridx(1);
      myc[2*dolog+iF]->SetGridy(1);
      //if (dolog==0) 
      myc[2*dolog+iF]->SetLogy(1);
      hP->SetMarkerStyle(20);
      hP->SetMarkerColor(2);
      hP->SetLineColor(2);
      hP->SetFillColor(2);
      hP->SetFillStyle(3004);
      hN->SetLineColor(4);
      hN->SetLineWidth(2);
      hN->SetMarkerStyle(21);
      hN->SetMarkerColor(4);
      if (dolog==0){
	hP->Rebin(5);
	hPi->Rebin(5);
	hN->Rebin(5);
      }
      //double max = 1.1*std::max(hP->GetMaximum(),hN->GetMaximum());
      //double max = 1.1*std::max(hPi->GetMaximum(),std::max(hP->GetMaximum(),hN->GetMaximum()));
      //hP->SetMaximum(max);
      if (dolog) hP->GetYaxis()->SetRangeUser(0.9,3.*1e6);
      else hP->GetYaxis()->SetRangeUser(0.9,3.*1e7);
      if (dolog==0) hP->GetXaxis()->SetRangeUser(0,30);
      hP->GetYaxis()->SetTitle("Particles");
      hP->Draw("hist");
      hP->Draw("PEsame");
      hN->Draw("histsame");
      hN->Draw("PEsame");
      
      if (iF==1){
	//hSi->SetMarkerStyle(23);
	//hSi->SetMarkerColor(5);
      hSi->SetLineColor(5);
      hSi->SetFillColor(5);
      hSi->SetFillStyle(3004);
      hSi->Draw("histsame");
      //hSi->Draw("PEsame");
      }

      hPi->SetMarkerStyle(22);
      hPi->SetMarkerColor(3);
      hPi->SetLineColor(3);
      //hPi->SetFillColor(3);
      //hPi->SetFillStyle(3005);
      hPi->Draw("histsame");
      hPi->Draw("PEsame");

      TLegend *leg = new TLegend(0.56,0.8,1.,1.);
      //if (dolog==1 && iF==1) leg = new TLegend(0.15,0.7,0.6,0.9);
      //else leg = new TLegend(0.56,0.8,1.,1.);
      leg->SetFillColor(10);
      std::ostringstream leglabel;
      leglabel << "n(protons) = " << hP->GetEntries();
      leg->AddEntry(hP,leglabel.str().c_str(),"PF");
      leglabel.str("");
      leglabel << "n(neutrons) = " << hN->GetEntries();
      leg->AddEntry(hN,leglabel.str().c_str(),"PL");
      leglabel.str("");
      leglabel << "n(#pi^{#pm}) = " << hPi->GetEntries();
      leg->AddEntry(hPi,leglabel.str().c_str(),"PL");
      if (iF==1){
	leglabel.str("");
	leglabel << "n(ions) = " << hSi->GetEntries();
	leg->AddEntry(hSi,leglabel.str().c_str(),"F");
      }
      leg->Draw("same");
      
      lat.SetTextSize(0.04);
      lat.DrawLatexNDC(0.1,0.96,"Pythia Minbias (30k)");
      if (iF==1) {
	lat.SetTextSize(0.03);
	lat.DrawLatexNDC(0.02,0.02,"HGCAL Standalone simulation");
	lat.SetTextSize(0.04);
      }
      
      //if (dolog==0){
	if (iF==0) {
	  lat.DrawLatexNDC(0.83,0.75,"@IP");
	  lat.DrawLatexNDC(0.8,0.69,"|#eta|<2.5");
	} else {
	  lat.DrawLatexNDC(0.52,0.75,"@shower max (layer 30)");
	  lat.DrawLatexNDC(0.73,0.69,"2.8 < #eta <3.0");
	  char buf[100];
	  sprintf(buf,"Si layer %d",siLayer);
	  lat.DrawLatexNDC(0.73,0.63,buf);
	}
	//} else {
	//if (iF==0) {
	//lat.DrawLatexNDC(0.15,0.45,"@IP");
	//lat.DrawLatexNDC(0.15,0.4,"|#eta|<2.5");
	//} else {
	//lat.DrawLatexNDC(0.15,0.45,"@shower max (layer 30)");
	//lat.DrawLatexNDC(0.15,0.4,"2.8 < #eta <3.0");
	//}
	//}

      myc[2*dolog+iF]->Update();
      myc[2*dolog+iF]->Print(("PLOTS/"+label[2*dolog+iF]+silabel+".pdf").c_str());
      myc[2*dolog+iF]->Print(("PLOTS/"+label[2*dolog+iF]+silabel+".png").c_str());

    }//dolog
  }//loop on files

  myc[4]->cd();
  gPad->SetLogy(1);
  gr->GetYaxis()->SetTitle("<particles>/event");
  gr->SetMarkerStyle(20);
  gr->SetFillColor(4);
  gr->Draw("AB1");
  std::vector<std::string> part;
  part.push_back("#bar{#Xi^{0}}");
  part.push_back("#Xi^{+}");
  part.push_back("#bar{#Sigma^{+}}");
  part.push_back("#bar{#Lambda}");
  part.push_back("#bar{#Sigma^{-}}");
  part.push_back("#bar{p}");
  part.push_back("#bar{n}");
  part.push_back("K^{-}");
  part.push_back("#pi^{-}");
  part.push_back("#bar{#nu_{#tau}}");
  part.push_back("#bar{#nu_{#mu}}");
  part.push_back("#mu^{+}");
  part.push_back("#bar{#nu_{e}}");
  part.push_back("e^{+}");
  part.push_back("e^{-}");
  part.push_back("#nu_{e}");
  part.push_back("#mu^{-}");
  part.push_back("#nu_{#mu}");
  part.push_back("#nu_{#tau}");
  part.push_back("#gamma");
  part.push_back("#pi^{0}");
  if (siLayer==1 || siLayer==3) part.push_back("#rho^{0}");
  part.push_back("K^{0}_{L}");
  part.push_back("#pi^{+}");
  part.push_back("#eta");
  part.push_back("K^{0}_{S}");
  part.push_back("K^{+}");
  if (siLayer!=1) part.push_back("#eta\'");
  part.push_back("n");
  part.push_back("p");
  part.push_back("#Sigma^{-}");
  part.push_back("#Lambda");
  part.push_back("#Sigma^{+}");
  part.push_back("#Xi^{-}");
  part.push_back("#Xi^{0}");
  part.push_back("Ions");

  TAxis *ax = gr->GetHistogram()->GetXaxis();
  double x1 = ax->GetBinLowEdge(1);
  double x2 = gr->GetN();
  gr->GetHistogram()->GetXaxis()->Set(gr->GetN(),x1,x2);
      
  std::cout << " parts n=" << part.size() << " axis n = " << ax->GetNbins() << " gr n=" << gr->GetN() 
	    << " x1=" << x1 << " x2=" << x2 
	    << std::endl;

  for(unsigned k=0;k<gr->GetN();k++){
    if (k<part.size()) gr->GetHistogram()->GetXaxis()->SetBinLabel(k+1,part[k].c_str());
  }
  lat.SetTextSize(0.04);
  lat.DrawLatexNDC(0.1,0.96,"Pythia Minbias");
  lat.SetTextSize(0.03);
  lat.DrawLatexNDC(0.02,0.02,"HGCAL Standalone simulation");
  lat.SetTextSize(0.04);
  lat.DrawLatexNDC(0.16,0.85,"@shower max (layer 30)");
  lat.DrawLatexNDC(0.16,0.78,"2.8 < #eta <3.0");
  char buf[100];
  sprintf(buf,"Si layer %d",siLayer);
  lat.DrawLatexNDC(0.16,0.70,buf);

  myc[4]->Update();
  myc[4]->Print(("PLOTS/ShowerMaxParticles"+silabel+".pdf").c_str());
  myc[4]->Print(("PLOTS/ShowerMaxParticles"+silabel+".png").c_str());


  myc[5]->cd();
  gPad->SetLogz(1);
  gPad->SetRightMargin(0.13);
  //hTvsLogp->GetYaxis()->SetTitle("");
  hTvsLogp->GetZaxis()->SetTitleOffset(0.75);
  hTvsLogp->Draw("colz");

  myc[5]->Update();
  myc[5]->Print(("PLOTS/ShowerMaxNeutronTimevsLogp"+silabel+".pdf").c_str());
  myc[5]->Print(("PLOTS/ShowerMaxNeutronTimevsLogp"+silabel+".png").c_str());

  return 0;

}//main
