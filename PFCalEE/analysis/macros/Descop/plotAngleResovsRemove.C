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

#include "../TDRStyle.h"

int plotAngleResovsRemove(){

  const unsigned nRemove = 13;
  //unsigned lToRemove[nRemove] = {27,24,21,8,19,3,6,1,10,14,17,12};

  TGraphErrors *grRes = new TGraphErrors();
  TGraphErrors *grResCheckF = new TGraphErrors();
  TGraphErrors *grResCheckB = new TGraphErrors();
  TGraphErrors *grResTP = new TGraphErrors();


  const unsigned nCanvas = 7;//nRemove;  
  TCanvas *mycR[nCanvas];
  for (unsigned iC(0);iC<nCanvas;++iC){
    std::ostringstream lName;
    lName.str("");
    lName << "mycR" << iC;
    mycR[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1500,1000);
    if (iC==1) mycR[iC]->Divide(4,4);
    if (iC>=2) mycR[iC]->Divide(7,2);
  }

  const unsigned nF = 6;
  TFile *fin[nF];
  std::ostringstream lname;
  std::ostringstream save;
  lname.str("");
  lname << "/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalDescopPos/gitV02-01-01/version12/gamma/v5_30_/eta16_et70_pu0_IC3.root";
  fin[0] = TFile::Open(lname.str().c_str());
  lname.str("");
  lname << "/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalDescopPosCheckFront/gitV02-01-01/version12/gamma/v5_30_/eta16_et70_pu0_IC3.root";
  fin[1] = TFile::Open(lname.str().c_str());
  lname.str("");
  lname << "/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalDescopPosCheckBack/gitV02-01-01/version12/gamma/v5_30_/eta16_et70_pu0_IC3.root";
  fin[2] = TFile::Open(lname.str().c_str());
  lname.str("");
  lname << "/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalDescopPos/gitV05-02-04/version30/gamma/fixedVtx_/eta17_et70_pu0_IC3.root";
  fin[3] = TFile::Open(lname.str().c_str());
  lname.str("");
  lname << "/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalDescopPos/gitV05-02-04/version34/gamma/fixedVtx_/eta17_et70_pu0_IC3.root";
  fin[4] = TFile::Open(lname.str().c_str());
  lname.str("");
  lname << "/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalDescopPos/gitV05-02-04/version35/gamma/fixedVtx_/eta17_et70_pu0_IC3.root";
  fin[5] = TFile::Open(lname.str().c_str());

  double val0 = 0;

  for (unsigned iF(0);iF<nF;++iF){
    fin[iF]->cd("PositionFit");
    TTree *tree = (TTree*)gDirectory->Get("PosFit");
    if (!tree) return 1;
    for (unsigned r(0); r<nRemove;++r){//loop on nremove
      if (iF>2 && r!=0) continue;
      if (iF==0) mycR[1]->cd(r+1);
      else if (iF==1) mycR[2]->cd(r+1);
      else if (iF==2) mycR[3]->cd(r+1);
      else mycR[1]->cd(iF+12);
      TH1F *hangle = new TH1F("hangle","",100,-0.1,0.1);
      lname.str("");
      lname << "recoEta-truthEta>>hangle";
      std::ostringstream lCut;
      lCut << "nRemove==" << r;
      tree->Draw(lname.str().c_str(),lCut.str().c_str()); 
      hangle->Fit("gaus");
      TF1 *fit = (TF1*)hangle->GetFunction("gaus");

      std::cout << " File " << iF << " r " << r << " " << fit->GetParameter(2) << "+/-" << fit->GetParError(2) << std::endl;

      if (iF==0){
	if (r==0) val0 = fit->GetParameter(2);
	grRes->SetPoint(r,r,fit->GetParameter(2));
	grRes->SetPointError(r,0,fit->GetParError(2));
      } else if (iF==1){
	grResCheckF->SetPoint(r,r,fit->GetParameter(2));
	grResCheckF->SetPointError(r,0,fit->GetParError(2));
      } else if (iF==2){
	grResCheckB->SetPoint(r,r,fit->GetParameter(2));
	grResCheckB->SetPointError(r,0,fit->GetParError(2));
      } else {
	double xval=2;
	if (iF==4) xval = 6;
	else if (iF==5) xval = 12;
	grResTP->SetPoint(iF-3,xval,fit->GetParameter(2));
	grResTP->SetPointError(iF-3,0,fit->GetParError(2));
      }

      if (iF<3){
	if (iF==0) mycR[4]->cd(r+1);
	else if (iF==1)  mycR[5]->cd(r+1);
	else if (iF==2)  mycR[6]->cd(r+1);
	tree->Draw("nLayersFit",lCut.str().c_str()); 
      }

    } 
   }//loop on files

  mycR[0]->cd();
  grRes->SetLineColor(1);
  grRes->SetMarkerColor(1);
  grRes->SetMarkerStyle(20);
  grRes->SetMinimum(0);
  grRes->SetMaximum(0.022);
  grRes->SetTitle(";n_{removed};#sigma_{#eta}");
  grRes->Draw("APEL");
  
  grResTP->SetLineColor(6);
  grResTP->SetMarkerColor(6);
  grResTP->SetMarkerStyle(21);
  grResTP->Draw("PEL");

  grResCheckF->SetLineColor(4);
  grResCheckF->SetMarkerColor(4);
  grResCheckF->SetMarkerStyle(22);
  grResCheckF->Draw("PEL");

  grResCheckB->SetLineColor(3);
  grResCheckB->SetMarkerColor(3);
  grResCheckB->SetMarkerStyle(23);
  grResCheckB->Draw("PEL");

  TF1 *ref = new TF1("ref","[0]*sqrt(30./(30-x))",0,13);
  ref->SetParameter(0,val0);
  ref->SetLineColor(7);
  ref->Draw("same");
  TLegend *leg = new TLegend(0.7,0.15,0.89,0.55);
  leg->SetFillColor(10);
  leg->AddEntry(grRes,"CMSSW","P");
  leg->AddEntry(grResCheckF,"CMSSW-drop-front","P");
  leg->AddEntry(grResCheckB,"CMSSW-drop-back","P");
  leg->AddEntry(grResTP,"TP","P");
  leg->AddEntry(ref,"#sqrt{#frac{30}{30-x}}","L");
  leg->Draw("same");

  mycR[0]->Update();
  save.str("");
  save << "PLOTS/RemoveScenarios_AngularResolution.pdf";
  mycR[0]->Print(save.str().c_str());
  
  mycR[1]->Update();
  save.str("");
  save << "PLOTS/RemoveScenarios_AngularResolutionFits.pdf";
  mycR[1]->Print(save.str().c_str());
  
  return 0;
  
}//main
