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

#include "TDRStyle.h"

int plotResoVsIC(){

  SetTdrStyle();


  const unsigned nIC = 10;
  const unsigned ICval[nIC] = {0,1,2,3,4,5,10,15,20,50};

  bool plotSR7 = true;

  std::ostringstream label;
  TFile *fcalib[nIC][2];
  
  TGraphErrors *constant = new TGraphErrors();
  constant->SetName("constant");
  constant->SetTitle(";intercalib. smearing");
  constant->SetMarkerStyle(20);
  constant->SetMarkerColor(1);
  constant->SetLineColor(1);
  TGraphErrors *constantSR7 =  new TGraphErrors();
  constantSR7->SetName("constantSR7");
  constantSR7->SetTitle(";intercalib. smearing");
  constantSR7->SetMarkerStyle(23);
  constantSR7->SetMarkerColor(2);
  constantSR7->SetLineColor(2);

  TGraphErrors *noise = (TGraphErrors *) constant->Clone("noise");
  TGraphErrors *sampling = (TGraphErrors *) constant->Clone("sampling");
  TGraphErrors *samplingSR7 = (TGraphErrors *) constantSR7->Clone("samplingSR7");

  TCanvas *mycReso = new TCanvas("mycReso","mycReso",1500,1000);
  mycReso->Divide(2,5);
  TCanvas *mycResoSR7 = new TCanvas("mycResoSR7","mycResoSR7",1500,1000);
  mycResoSR7->Divide(2,5);
  TCanvas *mycCalib = new TCanvas("mycCalib","mycCalib",1500,1000);
  mycCalib->Divide(2,5);
  TCanvas *mycCalibSR7 = new TCanvas("mycCalibSR7","mycCalibSR7",1500,1000);
  mycCalibSR7->Divide(2,5);
  TCanvas *mycR = new TCanvas("mycR","Sampling",1500,1000);
  TCanvas *mycC = new TCanvas("mycC","Constant",1500,1000);
  TCanvas *mycN = new TCanvas("mycN","Noise",1500,1000);

  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(0);

  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.5);

  TLatex lat;
  char buf[500];

  TGraphErrors *gr[nIC][2];
  TGraphErrors *grCalib[nIC][2];
  double y0=0;
  double y0_7=0;

  for (unsigned ic(0);ic<nIC;++ic){//loop on intercalib
    label.str("");
    label << "PLOTS/CalibReso";
    label << "_vsE";
    label << "_IC" << ICval[ic];
    label << ".root";
    fcalib[ic][0] = TFile::Open(label.str().c_str());
    if (!fcalib[ic][0]) {
      std::cout << " -- failed to open file: " << label.str() << std::endl;
      continue;
    }
    else {
      std::cout << " -- file " << label.str() << " successfully opened." << std::endl;
    }
    label.str("");
    label << "PLOTS/CalibReso";
    label << "_vsE";
    label << "_SR7_IC" << ICval[ic];
    label << "_try0";
    label << ".root";
    fcalib[ic][1] = TFile::Open(label.str().c_str());
    if (!fcalib[ic][1]) {
      std::cout << " -- failed to open file: " << label.str() << std::endl;
      continue;
    }
    else {
      std::cout << " -- file " << label.str() << " successfully opened." << std::endl;
    }
    fcalib[ic][0]->cd("SR2");
    gr[ic][0] = (TGraphErrors *)gDirectory->Get("resoRecoFit2eta21pu1");
    grCalib[ic][0] = (TGraphErrors *)gDirectory->Get("calibRecoDelta2eta21pu1");
    fcalib[ic][1]->cd("SR7");
    gr[ic][1] = (TGraphErrors *)gDirectory->Get("resoRecoFit7eta21pu1");
    grCalib[ic][1] = (TGraphErrors *)gDirectory->Get("calibRecoDelta7eta21pu1");


    TF1 *fit = gr[ic][0]->GetFunction("reso");
    TF1 *fit7 = gr[ic][1]->GetFunction("reso");
    mycReso->cd(ic+1);
    gr[ic][0]->Draw("APE");
    fit->SetLineColor(6);
    fit->Draw("same");
    lat.SetTextSize(0.1);
    sprintf(buf,"Single #gamma, #eta=2.1, 3#times3 cm^{2}");
    lat.DrawLatexNDC(0.2,0.8,buf);
    sprintf(buf,"ICsmear = %d %%",ICval[ic]);
    lat.DrawLatexNDC(0.2,0.7,buf);

    mycResoSR7->cd(ic+1);
    gr[ic][1]->Draw("APE");
    fit7->SetLineColor(6);
    fit7->Draw("same");
    lat.SetTextSize(0.1);
    sprintf(buf,"Single #gamma, #eta=2.1, all cells");
    lat.DrawLatexNDC(0.2,0.8,buf);
    sprintf(buf,"ICsmear = %d %%",ICval[ic]);
    lat.DrawLatexNDC(0.2,0.7,buf);


    mycCalib->cd(ic+1);
    gPad->SetGridy(1);
    grCalib[ic][0]->Draw("APE");
    lat.SetTextSize(0.1);
    sprintf(buf,"Single #gamma, #eta=2.1, 3#times3 cm^{2}");
    lat.DrawLatexNDC(0.2,0.8,buf);
    sprintf(buf,"ICsmear = %d %%",ICval[ic]);
    lat.DrawLatexNDC(0.2,0.7,buf);

    mycCalibSR7->cd(ic+1);
    gPad->SetGridy(1);
    grCalib[ic][1]->Draw("APE");
    lat.SetTextSize(0.1);
    sprintf(buf,"Single #gamma, #eta=2.1, all cells");
    lat.DrawLatexNDC(0.2,0.8,buf);
    sprintf(buf,"ICsmear = %d %%",ICval[ic]);
    lat.DrawLatexNDC(0.2,0.7,buf);

    if (ic==0){
      y0 = fit->GetParameter(1);
      y0_7 = fit7->GetParameter(1);
    }

    double cval = sqrt(pow(fit->GetParameter(1),2)-pow(y0,2));
    constant->SetPoint(ic,ICval[ic]/100.,cval);
    if (ic==0) constant->SetPointError(ic,0,fit->GetParError(1));
    else constant->SetPointError(ic,0,fit->GetParameter(1)*fit->GetParError(1)/cval);
    noise->SetPoint(ic,ICval[ic]/100.,fit->GetParameter(2));
    noise->SetPointError(ic,0,fit->GetParError(2));
    sampling->SetPoint(ic,ICval[ic]/100.,fit->GetParameter(0));
    sampling->SetPointError(ic,0,fit->GetParError(0));
    cval = sqrt(pow(fit7->GetParameter(1),2)-pow(y0_7,2));
    constantSR7->SetPoint(ic,ICval[ic]/100.,cval);
    if (ic==0) constantSR7->SetPointError(ic,0,fit7->GetParError(1));
    else constantSR7->SetPointError(ic,0,fit7->GetParameter(1)*fit7->GetParError(1)/cval);    
    //constantSR7->SetPoint(ic,ICval[ic]/100.,fit7->GetParameter(1));
    //constantSR7->SetPointError(ic,0,fit7->GetParError(1));
    samplingSR7->SetPoint(ic,ICval[ic]/100.,fit7->GetParameter(0));
    samplingSR7->SetPointError(ic,0,fit7->GetParError(0));

  }

  mycReso->Update();
  mycReso->Print("PLOTS/ResolutionFitvsIC_SR3.pdf");
  mycResoSR7->Update();
  mycResoSR7->Print("PLOTS/ResolutionFitvsIC_SR7.pdf");

  mycCalib->Update();
  mycCalib->Print("PLOTS/CalibrationRatiovsIC_SR3.pdf");
  mycCalibSR7->Update();
  mycCalibSR7->Print("PLOTS/CalibrationRatiovsIC_SR7.pdf");

  TLegend *leg = new TLegend(0.6,0.6,0.8,0.8);
  leg->SetFillColor(10);
  leg->AddEntry(sampling,"3#times3 cm^{2}","P");
  leg->AddEntry(samplingSR7,"All detector","P");
  mycR->cd();
  gPad->SetGridy(1);
  sampling->GetYaxis()->SetTitle("Sampling term (GeV^{1/2})");
  sampling->SetMinimum(0.2);
  sampling->SetMaximum(0.3);
  sampling->Draw("APE");
  samplingSR7->Draw("PEsame");
  lat.SetTextSize(0.04);
  sprintf(buf,"Single #gamma, #eta=2.1");
  lat.DrawLatexNDC(0.2,0.87,buf);
  lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");

  leg->Draw("same");
  mycR->Update();

  mycR->Print("PLOTS/SamplingvsIC.pdf");

  mycC->cd();
  gPad->SetLogx(1);
  gPad->SetGridy(1);
  gStyle->SetOptFit(0);
  //gStyle->SetStatH(0.1);
  //gStyle->SetStatW(0.2);

  constant->GetYaxis()->SetTitle("Constant from intercalib.");
  constant->SetMinimum(0);
  constant->SetMaximum(0.15);
  constant->Draw("APE");
  constantSR7->Draw("PEsame");
  sprintf(buf,"Single #gamma, #eta=2.1");
  lat.DrawLatexNDC(0.2,0.87,buf);
  lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");

  TF1 *BE = new TF1("BE","sqrt([0]*[0] + pow(x*1/sqrt([1]),2))",0,1);
  BE->SetParameters(0,30);
  //  BE->SetParLimits(0,1,1);
  BE->FixParameter(0,0);
  BE->SetLineColor(1);
  constant->Fit("BE","BI");

  lat.SetTextColor(1);
  //sprintf(buf,"c #propto c_{0} #oplus #frac{x}{#sqrt{n}}, n=%3.1f #pm %3.1f",BE->GetParameter(1),BE->GetParError(1));
  sprintf(buf,"c_{ic}=#frac{x}{#sqrt{n}}, n=%3.1f #pm %3.1f",BE->GetParameter(1),BE->GetParError(1));
  lat.DrawLatexNDC(0.2,0.77,"c=c_{0} #oplus c_{ic}");
  lat.DrawLatexNDC(0.2,0.67,buf);

  //BE->SetParameter(0,y0_7);
  BE->SetParameters(0,30);
  BE->FixParameter(0,0);
  BE->SetLineColor(2);
  constantSR7->Fit("BE","BI","same");
  //BE->Draw("same");
  lat.SetTextColor(2);
  //sprintf(buf,"c #propto c_{0} #oplus #frac{x}{#sqrt{n}}, n=%3.1f #pm %3.1f",BE->GetParameter(1),BE->GetParError(1));
  sprintf(buf,"c_{ic}=#frac{x}{#sqrt{n}}, n=%3.1f #pm %3.1f",BE->GetParameter(1),BE->GetParError(1));
  lat.DrawLatexNDC(0.2,0.6,buf);
  lat.SetTextColor(1);

  leg->Draw("same");
  mycC->Update();
  mycC->Print("PLOTS/ConstantvsIC.pdf");


  mycN->cd();
  noise->GetYaxis()->SetTitle("Noise term (GeV)");
  noise->Draw("APE");
  sprintf(buf,"Single #gamma, #eta=2.1, 3#times3 cm^{2}");
  lat.DrawLatexNDC(0.2,0.87,buf);
  lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");
  mycN->Print("PLOTS/NoisevsIC.pdf");

  




  return 0;

}//main
