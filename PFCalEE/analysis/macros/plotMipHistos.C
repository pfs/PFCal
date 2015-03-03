#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>

#include "TFile.h"
#include "TTree.h"
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
#include "TMath.h"
#include "TVectorD.h"
#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"
#include "TSystem.h"

#include "TDRStyle.h"


Double_t langaufun(Double_t *x, Double_t *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation), 
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = 0;//-0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;


      // MP shift correction
      mpc = par[1] - mpshift * par[0]; 

      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }

      return (par[2] * step * sum * invsq2pi / par[3]);
};


Double_t langaufungaus(Double_t *x, Double_t *par) {
  
  return langaufun(x,par)+par[4]*TMath::Gaus(x[0],par[5],par[6],true);

};

int plotMipHistos(){

  ///////////////////////////////////////////////////////////////////
  ////////// Hardcoded config ///////////////////////////////////////
  ///////////////////////////////////////////////////////////////////

  const unsigned nEta = 1;
  //const unsigned nNoise = 10;
  const unsigned nNoise = 6;//10;

  //const double noise[nNoise] = {0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5};
  const double noise[nNoise] = {0,0.6,0.7,0.8,0.9,1.0};

  const double eta[nEta] = {2.85};//1.7,2.0,2.5};

  const double deta = 0.05;

  const unsigned nLayers = 30;

  //bool doSignalOnly = true;
  bool doSignalOnly = false;

  bool oneOnly = true;
  const double Ethresh = 0.9;
  const double EthreshMax = 5.0;
  const double EmaxCut = 0.05;
  const unsigned layerRange = 2;

  std::ostringstream suffix;
  suffix << "Eta" << eta[0];
  suffix << "_thresh" << Ethresh ;
  suffix << "_" << EthreshMax;
  suffix << "_EmaxNeighbour" << EmaxCut;
  suffix << "_trk" << 1+2*layerRange << "layers";
  suffix << "_highnoise";
  //suffix << "_1x1";  
  if (oneOnly) suffix << "_onlyOne";

  ///////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////

  std::ostringstream outFilePath;
  outFilePath << "/afs/cern.ch/work/a/amagnan/PFCalEEAna/PLOTS/gitV00-02-12/version12/MinBias/Histos" << suffix.str() << ".root";
  TFile *file = TFile::Open(outFilePath.str().c_str());

  TFile *file0 = TFile::Open("/afs/cern.ch/work/a/amagnan/PFCalEEAna/PLOTS/gitV00-02-12/version12/MinBias/HistosEta2.85_thresh0.9_5_EmaxNeighbour0.05_trk5layers_onlyOne.root");

  if (doSignalOnly) suffix << "_sigOnly";


  TF1 *landaugaus = new TF1("landaugaus",langaufungaus,0,5,7);

  double mpshift  = 0;//0.22278298;

  TH1F *hitSpectrum[nEta][nNoise][nLayers];
  TH1F *hitSpectrumSig[nEta][nNoise][nLayers];
  TH1F *p_nHits[nEta][nNoise][nLayers];

  TH1F *p_EmaxNeighbour[nEta][nNoise];

  std::ostringstream label;

  SetTdrStyle();
  TCanvas *myc = new TCanvas("myc","myc",2000,1000);
  TCanvas *mycP = new TCanvas("mycP","mycP",1);
  TCanvas *mycS = new TCanvas("mycS","mycS",2000,1000);
  mycS->Divide(3,2);
  TCanvas *mycM = new TCanvas("mycM","mycM",2000,1000);
  mycM->Divide(3,2);
  TCanvas *mycH = new TCanvas("mycH","mycH",1500,1000);
  TCanvas *mycL = new TCanvas("mycL","mycL",2000,1000);
  mycL->Divide(3,2);


  // TCanvas *mycL[nLayers];
  // for (unsigned il(0); il<nLayers;++il){
  //   label.str("");
  //   label << "mycL" << il ;
  //   mycL[il] = new TCanvas(label.str().c_str(),label.str().c_str(),1500,1000);
  //   mycL[il]->Divide(3,2);
  // }


  gStyle->SetOptStat("eMRuo");
  gStyle->SetOptFit(1111);
  TLatex lat;
  char buf[500];

  //TH1F *hitSpectrum[nEta][nNoise][nLayers];
  TH1F *hitSpectrumAll[nEta][nNoise];
  TH1F *hitSpectrumAllSig[nEta][nNoise];

  TGraphErrors *grNhits[nEta][nNoise];
  TGraphErrors *grMPV[nEta][nNoise];
  TGraphErrors *grSigma[nEta][nNoise];

  TLegend *leg = new TLegend(0.7,0.7,0.94,0.94);
  leg->SetFillColor(0);
  TLegend *legN = new TLegend(0.8,0.55,0.94,0.94);
  legN->SetFillColor(0);

  for (unsigned ie(0);ie<nEta;++ie){

    for (unsigned in(0); in<nNoise;++in){
      grNhits[ie][in] = new TGraphErrors();
      label.str("");
      label << "grNhits_eta" << eta[ie]*10 << "_noise" << noise[in]*100;
      grNhits[ie][in]->SetName(label.str().c_str());
      grMPV[ie][in] = new TGraphErrors();
      label.str("");
      label << "grMPV_eta" << eta[ie]*10 << "_noise" << noise[in]*100;
      grMPV[ie][in]->SetName(label.str().c_str());
      grSigma[ie][in] = new TGraphErrors();
      label.str("");
      label << "grSigma_eta" << eta[ie]*10 << "_noise" << noise[in]*100;
      grSigma[ie][in]->SetName(label.str().c_str());

      if (in==0) file0->cd();
      else file->cd();
      label.str("");
      label << "EmaxNeighbour_eta" << eta[ie]*10 << "_noise" << noise[in]*100;
      p_EmaxNeighbour[ie][in] = (TH1F*)gDirectory->Get(label.str().c_str());
      for (unsigned il(0); il<nLayers;++il){
	label.str("");
	label << "hitSpectrum";
	if (doSignalOnly) label << "Sig";
	label << "_eta" << eta[ie]*10 << "_noise" << noise[in]*100 << "_layer" << il;
	hitSpectrum[ie][in][il] = (TH1F*)gDirectory->Get(label.str().c_str());
	if (il==0) hitSpectrumAll[ie][in] = (TH1F*)hitSpectrum[ie][in][il]->Clone(label.str().c_str());
	else hitSpectrumAll[ie][in]->Add(hitSpectrum[ie][in][il]);

	label.str("");
	label << "hitSpectrumSig_eta" << eta[ie]*10 << "_noise" << noise[in]*100 << "_layer" << il;
	hitSpectrumSig[ie][in][il] = (TH1F*)gDirectory->Get(label.str().c_str());
 	if (il==0) hitSpectrumAllSig[ie][in] = (TH1F*)hitSpectrumSig[ie][in][il]->Clone(label.str().c_str());
	else hitSpectrumAllSig[ie][in]->Add(hitSpectrumSig[ie][in][il]);

	label.str("");
	label << "Nhits_eta" << eta[ie]*10 << "_noise" << noise[in]*100 << "_layer" << il;
	p_nHits[ie][in][il] = (TH1F*)gDirectory->Get(label.str().c_str());
	
      }
    }
  }

  myc->cd();
  gStyle->SetOptStat(0);
  myc->SetLogy(1);
  for (unsigned ie(0);ie<nEta;++ie){

    for (unsigned in(0); in<nNoise;++in){
      p_EmaxNeighbour[ie][in]->SetLineColor(in+1);
      p_EmaxNeighbour[ie][in]->SetMarkerColor(in+1);
      p_EmaxNeighbour[ie][in]->SetMarkerStyle(ie+20);
      p_EmaxNeighbour[ie][in]->SetMaximum(100000);//p_EmaxNeighbour[0][nNoise-1]->GetMaximum());
      if (ie==0 && in==0) p_EmaxNeighbour[ie][in]->Draw("PEL");
      else p_EmaxNeighbour[ie][in]->Draw("PELsame");
      label.str("");
      label << "#sigma_{N} = " << noise[in] << " mips" ;
      if (ie==0) legN->AddEntry(p_EmaxNeighbour[ie][in],label.str().c_str(),"P");
    }
  }
  legN->Draw("same");
  myc->Update();
  label.str("");
  label << "PLOTS/" << suffix.str() << "/EmaxNeighbour";
  //label << suffix.str();
  label << ".pdf";
  myc->Print(label.str().c_str());
  //return 1;
  myc->Clear();
  myc->SetLogy(0);
  myc->Divide(3,2);

  gStyle->SetOptStat("eMRuo");

  for (unsigned ie(0);ie<nEta;++ie){
    for (unsigned il(0); il<nLayers;++il){
      for (unsigned in(0); in<nNoise;++in){
	mycL->cd(in+1);
	if (noise[in]>0.3) hitSpectrum[ie][in][il]->Rebin(2);
	hitSpectrum[ie][in][il]->Draw();
	if (hitSpectrum[ie][in][il]->GetEntries()==0) continue;
	TF1 *fit = 0;
	double mpc = 0;
	double minfit = 0.5;
	double maxfit = noise[in]<0.31 ? 1.5 : 2.5;//1.7;
	if (in==0){
	  hitSpectrum[ie][in][il]->Fit("landau","R+","same",minfit,maxfit);
	  fit = (TF1*)hitSpectrum[ie][in][il]->GetFunction("landau");
	  if (!fit) continue;
	  landaugaus->SetParameters(fit->GetParameter(2),fit->GetParameter(1),fit->GetParameter(0),noise[in],fit->GetParameter(0),0,noise[in]);
	  if (doSignalOnly) {
	    landaugaus->SetParameter(6,0);
	    landaugaus->SetParameter(4,0);
	    landaugaus->FixParameter(4,0);
	  }
	  mpc = fit->GetParameter(1) - mpshift * fit->GetParameter(2); 
	}
	else {
	  landaugaus->FixParameter(3,noise[in]);
	  landaugaus->SetParLimits(0,0.,0.2);
	  landaugaus->FixParameter(5,0);
	  if (doSignalOnly) landaugaus->FixParameter(6,0.);
	  else {
	    landaugaus->FixParameter(6,noise[in]);
	    minfit = 0.;
	  }
	  hitSpectrum[ie][in][il]->Fit("landaugaus","BR+","same",minfit,maxfit);
	  fit = (TF1*)hitSpectrum[ie][in][il]->GetFunction("landaugaus");
	  if (!fit) continue;
	  mpc = fit->GetParameter(1);
	}
	sprintf(buf,"Layer %d, Noise = %3.2f MIPs",il,noise[in]);
	lat.DrawLatexNDC(0.2,0.85,buf);
	sprintf(buf,"%3.2f < #eta < %3.2f",eta[ie]-deta,eta[ie]+deta);
	lat.DrawLatexNDC(0.2,0.75,buf);


	grNhits[ie][in]->SetPoint(il,il,p_nHits[ie][in][il]->GetMean());
	grNhits[ie][in]->SetPointError(il,0,p_nHits[ie][in][il]->GetRMS());

	grMPV[ie][in]->SetPoint(il,il,mpc);
	grMPV[ie][in]->SetPointError(il,0,fit->GetParError(1));
	if (in==0){
	  grSigma[ie][in]->SetPoint(il,il,fit->GetParameter(2));
	  grSigma[ie][in]->SetPointError(il,0,fit->GetParError(2));
	}
	else {
	  grSigma[ie][in]->SetPoint(il,il,fit->GetParameter(0));
	  grSigma[ie][in]->SetPointError(il,0,fit->GetParError(0));
	}
      }
      
      label.str("");
      label << "PLOTS/" << suffix.str() << "/MipSpectrum_eta" << eta[ie]*10 ;
      label << "_layer" << il;
      //label << suffix.str();
      label << ".pdf";
      mycL->Print(label.str().c_str());

    }

    for (unsigned in(0); in<nNoise;++in){
      myc->cd(in+1);
      gPad->SetLogy(1);
      if (noise[in]>0.3) {
	hitSpectrumAll[ie][in]->Rebin(2);
	hitSpectrumAllSig[ie][in]->Rebin(2);
      }
      hitSpectrumAll[ie][in]->Draw();
      hitSpectrumAllSig[ie][in]->SetLineColor(4);
      hitSpectrumAllSig[ie][in]->Draw("same");
      double minfit = 0.5;
      double maxfit = noise[in]<0.31 ? 1.5 : 2.5;//1.7;
      if (in==0){
	minfit = 0.75;
	maxfit = 1.4;
	hitSpectrumAll[ie][in]->Fit("landau","R+","same",minfit,maxfit);
	TF1 *fit = (TF1*)hitSpectrumAll[ie][in]->GetFunction("landau");
	if (!fit) continue;
	landaugaus->SetParameters(fit->GetParameter(2),fit->GetParameter(1),fit->GetParameter(0),noise[in],fit->GetParameter(0),0,noise[in]);
	if (doSignalOnly) {
	  landaugaus->SetParameter(6,0);
	  landaugaus->SetParameter(4,0);
	  landaugaus->FixParameter(4,0);
	}
      }
      else {
	//landaugaus->SetParameters(0.11,1,200000.,noise[in]);
	landaugaus->FixParameter(3,noise[in]);
	landaugaus->SetParLimits(0,0.,0.2);
	landaugaus->FixParameter(5,0);
	if (doSignalOnly) landaugaus->FixParameter(6,0);
	else {
	  landaugaus->FixParameter(6,noise[in]);
	  minfit = 0.;
	}
	//if (in>=nNoise-2) landaugaus->FixParameter(0,landaugaus->GetParameter(0));
	//landaugaus->SetLineColor(6);
	//landaugaus->Draw("same");
	hitSpectrumAll[ie][in]->Fit("landaugaus","BR+","same",minfit,maxfit);
      }
      sprintf(buf,"Noise = %3.2f MIPs",noise[in]);
      lat.DrawLatexNDC(0.2,0.85,buf);
      sprintf(buf,"%3.2f < #eta < %3.2f",eta[ie]-deta,eta[ie]+deta);
      lat.DrawLatexNDC(0.2,0.75,buf);

      lat.SetTextColor(4);
      lat.DrawLatexNDC(0.7,0.25,"Signal only");
      lat.SetTextColor(1);

      if (noise[in] == 0.3){
	mycP->cd();
	gStyle->SetOptFit(0);
	hitSpectrumAll[ie][in]->Draw();
	sprintf(buf,"Noise = %3.2f MIPs",noise[in]);
	lat.DrawLatexNDC(0.6,0.85,buf);
	sprintf(buf,"%3.2f < #eta < %3.2f",eta[ie]-deta,eta[ie]+deta);
	lat.DrawLatexNDC(0.6,0.75,buf);
	label.str("");
	label << "PLOTS/" << suffix.str() << "/MipSpectrum_eta" << eta[ie]*10 << "_noise" << noise[in] ;
	//label << suffix.str();
	//label << ".pdf";
	mycP->Print((label.str()+".pdf").c_str());
	mycP->Print((label.str()+".C").c_str());
      }
      gStyle->SetOptFit(1111);

    }
    label.str("");
    label << "PLOTS/" << suffix.str() << "/MipSpectrum_eta" << eta[ie]*10 ;
    //label << suffix.str();
    label << ".pdf";
    myc->Print(label.str().c_str());

    //return 1;

    for (unsigned in(0); in<nNoise;++in){
      myc->cd(in+1);
      gPad->SetLogy(0);
      hitSpectrumAll[ie][in]->Draw();
      hitSpectrumAllSig[ie][in]->Draw("same");
      sprintf(buf,"Noise = %3.1f MIPs",noise[in]);
      lat.DrawLatexNDC(0.35,0.8,buf);
      sprintf(buf,"%3.2f < #eta < %3.1f",eta[ie]-deta,eta[ie]+deta);
      lat.DrawLatexNDC(0.35,0.87,buf);
      lat.SetTextColor(4);
      lat.DrawLatexNDC(0.7,0.2,"From signal");
      lat.SetTextColor(1);
      lat.DrawLatexNDC(0.7,0.25,"All");
      lat.SetTextColor(2);
      lat.DrawLatexNDC(0.78,0.25,"+ fit");
      lat.SetTextColor(1);

      mycM->cd(in+1);
      gPad->SetGridy(1);
      grMPV[ie][in]->SetTitle(";layer;Landau MPV (MIPs)");
      grMPV[ie][in]->SetMarkerStyle(20+ie);
      grMPV[ie][in]->SetMarkerColor(1+ie);
      grMPV[ie][in]->SetLineColor(1+ie);
      grMPV[ie][in]->SetMinimum(0.8);
      grMPV[ie][in]->SetMaximum(1.2);
      grMPV[ie][in]->Draw(ie==0?"APL":"PLsame");
      label.str("");
      label << "#eta = " << eta[ie] ;
      if (in==0) leg->AddEntry(grMPV[ie][in],label.str().c_str(),"P");
      if (ie==nEta-1) leg->Draw("same");
      sprintf(buf,"Noise = %3.2f MIPs",noise[in]);
      lat.DrawLatexNDC(0.2,0.85,buf);

      mycS->cd(in+1);
      gPad->SetGridy(1);
      grSigma[ie][in]->SetTitle(";layer;Landau #sigma (MIPs)");
      grSigma[ie][in]->SetMarkerStyle(20+ie);
      grSigma[ie][in]->SetMarkerColor(1+ie);
      grSigma[ie][in]->SetLineColor(1+ie);
      grSigma[ie][in]->SetMinimum(0.);
      grSigma[ie][in]->SetMaximum(0.3);
      grSigma[ie][in]->Draw(ie==0?"APL":"PLsame");
      if (ie==nEta-1) leg->Draw("same");
      sprintf(buf,"Noise = %3.2f MIPs",noise[in]);
      lat.DrawLatexNDC(0.2,0.85,buf);

    }
    label.str("");
    label << "PLOTS/" << suffix.str() << "/MipSpectrum_eta" << eta[ie]*10 ;
    //label << suffix.str();
    label << "_nolog.pdf";
    myc->Print(label.str().c_str());

  }
  label.str("");
  label << "PLOTS/" << suffix.str() << "/LandauMPV";
  //label << suffix.str();
  label << ".pdf";
  mycM->Print(label.str().c_str());
  label.str("");
  label << "PLOTS/" << suffix.str() << "/LandauSigma";
  //label << suffix.str();
  label << ".pdf";
  mycS->Print(label.str().c_str());

  mycS->Clear();
  mycM->Clear();
  for (unsigned ie(0);ie<nEta;++ie){
    for (unsigned in(0); in<nNoise;++in){
      mycH->cd();
      gPad->SetGridy(1);
      grNhits[ie][in]->SetMarkerStyle(20+in);
      grNhits[ie][in]->SetMarkerColor(1+in);
      grNhits[ie][in]->SetLineColor(1+in);
      if (in==4) {
	grNhits[ie][in]->SetMarkerColor(6);
	grNhits[ie][in]->SetLineColor(6);
      }
      grNhits[ie][in]->SetMinimum(0);
      grNhits[ie][in]->SetMaximum(1000);
      grNhits[ie][in]->Draw(in==0?"APL":"PLsame");
      mycM->cd();
      gPad->SetGridy(1);
      grMPV[ie][in]->SetMarkerStyle(20+in);
      grMPV[ie][in]->SetMarkerColor(1+in);
      grMPV[ie][in]->SetLineColor(1+in);
      if (in==4) {
	grMPV[ie][in]->SetMarkerColor(6);
	grMPV[ie][in]->SetLineColor(6);
      }
      grMPV[ie][in]->SetMinimum(0.8);
      grMPV[ie][in]->SetMaximum(1.2);
      grMPV[ie][in]->Draw(in==0?"APL":"PLsame");
      //label.str("");
      //label << "#sigma_{N} = " << noise[in] << " mips" ;
      //if (ie==0) legN->AddEntry(grMPV[ie][in],label.str().c_str(),"P");
      mycS->cd();
      gPad->SetGridy(1);
      grSigma[ie][in]->SetMarkerStyle(20+in);
      grSigma[ie][in]->SetMarkerColor(1+in);
      grSigma[ie][in]->SetLineColor(1+in);
      if (in==4) {
	grSigma[ie][in]->SetMarkerColor(6);
	grSigma[ie][in]->SetLineColor(6);
      }
      grSigma[ie][in]->SetMinimum(0.);
      grSigma[ie][in]->SetMaximum(0.3);
      grSigma[ie][in]->Draw(in==0?"APL":"PLsame");

    }
    mycH->cd();
    legN->Draw("same");
    sprintf(buf,"%3.2f < #eta < %3.2f",eta[ie]-deta,eta[ie]+deta);
    lat.DrawLatexNDC(0.2,0.85,buf);
    label.str("");
    label << "PLOTS/" << suffix.str() << "/Nhits_eta" << eta[ie]*10 ;
    //label << suffix.str();
    label << ".pdf";
    mycH->Print(label.str().c_str());

    mycM->cd();
    legN->Draw("same");
    sprintf(buf,"%3.2f < #eta < %3.2f",eta[ie]-deta,eta[ie]+deta);
    lat.DrawLatexNDC(0.2,0.85,buf);
    label.str("");
    label << "PLOTS/" << suffix.str() << "/LandauMPV_eta" << eta[ie]*10 ;
    //label << suffix.str();
    label << ".pdf";
    mycM->Print(label.str().c_str());
    mycS->cd();
    legN->Draw("same");
    sprintf(buf,"%3.2f < #eta < %3.2f",eta[ie]-deta,eta[ie]+deta);
    lat.DrawLatexNDC(0.2,0.85,buf);
    label.str("");
    label << "PLOTS/" << suffix.str() << "/LandauSigma_eta" << eta[ie]*10 ;
    //label << suffix.str();
    label << ".pdf";
    mycS->Print(label.str().c_str());

  }

  return 0;
}//main
