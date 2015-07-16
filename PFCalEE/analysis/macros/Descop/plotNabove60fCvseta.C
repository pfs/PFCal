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

//#define ISPI

int plotNabove60fCvseta(){//main

  SetTdrStyle();
#ifdef ISPI
  bool isPi = true;
  std::string particle = "pi-";
#else
  bool isPi = false;
  std::string particle = "gamma";
#endif

#ifdef ISPI
  const unsigned npt = 15;
  const double pt[npt] = {1,2,3,4,5,6,7,8,9,10,12,14,16,18,20};
#else
  const unsigned npt = 9;
  const double pt[npt] = {3,5,10,20,30,50,70,100,200};
#endif

#ifdef ISPI

  const unsigned neta = 5;
  const double eta[neta] = {1.75,2.0,2.25,2.5,2.75};
#else
  const unsigned neta = 7;
  const double eta[neta] = {1.7,1.9,2.1,2.3,2.5,2.5,2.7};
#endif

  const std::string path = "/afs/cern.ch/work/a/amagnan/PFCalEEAna//HGCalTime/gitV05-02-04/"+particle+"/";

  const unsigned limit = 1;

  const unsigned nV = 3;
#ifdef ISPI
  const unsigned v[nV] = {33,36,37};
#else
  unsigned v[nV] = {30,34,35};
#endif

  TCanvas *myc2 = new TCanvas("myc2","myc2",1500,1000);
  TCanvas *myc[nV];
  TCanvas *myc1[nV];

  for (unsigned iv(0); iv<nV;++iv){//loop on version
    std::ostringstream lname;
    lname << "myc_" << v[iv]; 
    myc[iv] = new TCanvas(lname.str().c_str(),lname.str().c_str(),1500,1000);
    lname.str("");
    lname << "myc1_" << v[iv]; 
    myc1[iv] = new TCanvas(lname.str().c_str(),lname.str().c_str(),1500,1000);
  }
  TGraphErrors *grDummy = new TGraphErrors();
  grDummy->SetMinimum(0);
  grDummy->SetMaximum(1.1);
  TGraphErrors *grDummy1 = new TGraphErrors();
  grDummy1->SetMinimum(isPi ?0 : 0.5);
  grDummy1->SetMaximum(isPi ? 26 : 1000);
  for (unsigned ieta(0);ieta<neta;++ieta){//loop on pt
    //if (ieta==4) grDummy->SetPoint(ieta,eta[ieta]-0.05,0);
    //else if (ieta==5) grDummy->SetPoint(ieta,eta[ieta]+0.05,0);
    grDummy->SetPoint(ieta,eta[ieta],0);
    //if (ieta==4) grDummy1->SetPoint(ieta,eta[ieta]-0.05,0);
    //else if (ieta==5) grDummy1->SetPoint(ieta,eta[ieta]+0.05,0);
    grDummy1->SetPoint(ieta,eta[ieta],0);
  }
  std::ostringstream label;
  label << ";#eta;proba(n_{E>60fC}#geq" << limit << ")";
  grDummy->SetTitle(label.str().c_str());
  label.str("");
  label << ";#eta;<n_{E>60fC}>";
  grDummy1->SetTitle(label.str().c_str());

  TGraphErrors *gr[nV][npt];
  TGraphErrors *grMean[nV][npt];
  TTree *tree = 0;
  TFile *ftmp = 0;

  TLegend *leg;
  if (isPi) leg = new TLegend(0.72,0.43,0.87,0.93);
  else leg = new TLegend(0.78,0.11,0.94,0.54);
  leg->SetFillColor(10);
  for (unsigned iv(0); iv<nV;++iv){//loop on version
    myc[iv]->cd();
    grDummy->Draw("AL");
    myc1[iv]->cd();
    if (!isPi) gPad->SetLogy(1);
    gPad->SetGridy(1);
    grDummy1->Draw("AL");
    grDummy1->GetYaxis()->SetNdivisions(121);
    //if (!isPi) grDummy1->GetXaxis()->SetRangeUser(2.4,2.6);
 
    for (unsigned ipt(0); ipt<npt;++ipt){//loop on pt
      std::ostringstream lname;
      lname << "grv" << v[iv]<< "_pt" << pt[ipt];
      gr[iv][ipt] = new TGraphErrors();
      gr[iv][ipt]->SetName(lname.str().c_str());
      grMean[iv][ipt] = new TGraphErrors();
      lname << "_mean";
      grMean[iv][ipt]->SetName(lname.str().c_str());
      //set dummy values to have full range for sure...
      for (unsigned ieta(0);ieta<neta;++ieta){//loop on pt
	std::ostringstream fname;
	fname << path << "v" << v[iv] << "_et" << pt[ipt] << "_eta" << eta[ieta];
	if (eta[ieta]==2) fname << ".0";
	if (ieta==4) fname << "/NAbove26_old.root";
	else fname << "/NAbove26.root";
	ftmp = TFile::Open(fname.str().c_str());
	if (!ftmp) {
	  continue;
	}
	ftmp->cd();
	tree = (TTree*)gDirectory->Get("outtree");
	if (!tree) {
	  continue;
	}
	TH1F *hprob = new TH1F("hprob",";n(E>26 mips);proba",30,0,30);
	hprob->StatOverflows(1);
	hprob->Sumw2();
	myc2->cd();
	tree->Draw("nover>>hprob");
	unsigned nentries = hprob->GetEntries();
	std::cout << "v " << v[iv] 
		  << " et " << pt[ipt] 
		  << " eta " << eta[ieta];
	double neweta = eta[ieta];
	//if (ieta==4) neweta = 2.45;
	//else if (ieta==5) neweta = 2.55;
	if (nentries>0) {
	  std::cout << " mean " << hprob->GetMean() << " " << hprob->GetMeanError() ;
	  grMean[iv][ipt]->SetPoint(grMean[iv][ipt]->GetN(),neweta,hprob->GetMean());
	  grMean[iv][ipt]->SetPointError(grMean[iv][ipt]->GetN()-1,0,hprob->GetMeanError());
	  hprob->Scale(1./nentries);
	}
	else {
	  continue;
	}
	double error = 0;
	double integral = hprob->IntegralAndError(limit+1,31,error);
 
	std::cout << " entries=" << hprob->GetEntries() 
		  << " integral(" << limit+1 << ",N) " << integral << "+/-" << error
		  << std::endl;
	gr[iv][ipt]->SetPoint(gr[iv][ipt]->GetN(),neweta,integral);
	gr[iv][ipt]->SetPointError(gr[iv][ipt]->GetN()-1,0,error);
	hprob->Delete();
	ftmp->Close();
      }//loop on eta
      myc[iv]->cd();
      std::cout << "v " << v[iv] 
		<< " pt " << pt[ipt]
		<< " npoints = " << gr[iv][ipt]->GetN()
		<< " nmean = " << grMean[iv][ipt]->GetN()
		<< std::endl;

      gr[iv][ipt]->SetMarkerColor(ipt%8+(ipt%8<4?1:2));
      gr[iv][ipt]->SetLineColor(ipt%8+(ipt%8<4?1:2));
      gr[iv][ipt]->SetMarkerStyle(20+ipt);
      gr[iv][ipt]->SetMinimum(0.0);
      gr[iv][ipt]->SetMaximum(1.1);
      std::ostringstream label;
      label << ";#eta;proba(n_{E>60fC}#geq" << limit << ")";
      gr[iv][ipt]->SetTitle(label.str().c_str());
      gr[iv][ipt]->Draw("PE");
      label.str("");
      label << "p_{T}=" << pt[ipt] << " GeV";
      if (iv==0) leg->AddEntry(gr[iv][ipt],label.str().c_str(),"P");

      myc1[iv]->cd();
      grMean[iv][ipt]->SetMarkerColor(ipt%8+(ipt%8<4?1:2));
      grMean[iv][ipt]->SetLineColor(ipt%8+(ipt%8<4?1:2));
      grMean[iv][ipt]->SetMarkerStyle(20+ipt);
      grMean[iv][ipt]->SetMinimum(0.01);
      grMean[iv][ipt]->SetMaximum(25);
      label.str("");
      label << ";#eta;<n_{E>60fC}>";
      grMean[iv][ipt]->SetTitle(label.str().c_str());
      grMean[iv][ipt]->Draw("PEL");

    }//loop on pt
    myc[iv]->cd();
    leg->Draw("same");
    myc[iv]->Update();
    std::ostringstream lsave;
    lsave << "SummaryNabove60fC_" << particle << "_limit" << limit << "_vseta_version" << v[iv] << ".pdf";
    myc[iv]->Print(lsave.str().c_str());

    myc1[iv]->cd();
    leg->Draw("same");
    myc1[iv]->Update();
    lsave.str("");
    lsave << "SummaryNMeanabove60fC_" << particle << "_vseta_version" << v[iv] << ".pdf";
    myc1[iv]->Print(lsave.str().c_str());


  }//loop on version

  return 0;
}//main
