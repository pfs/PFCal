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

#define ISPI

int plotNabove60fCvspt(){//main

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
  const unsigned npt = 8;
  const double pt[npt] = {3,5,10,30,50,70,100,200};
#endif

#ifdef ISPI
  const unsigned neta = 5;
  const double eta[neta] = {1.75,2.0,2.25,2.5,2.75};
#else
  const unsigned neta = 1;
  const double eta[neta] = {2.5};
#endif

  const std::string path = "/afs/cern.ch/work/a/amagnan/PFCalEEAna//HGCalTime/gitV05-02-04/"+particle+"/";

  const unsigned limit = 1;

  const unsigned nV = 3;
#ifdef ISPI
  unsigned v[nV] = {33,36,37};
#else
  unsigned v[nV] = {30,34,35};
#endif

  TCanvas *myc2 = new TCanvas("myc2","myc2",1500,1000);
  TCanvas *myc[neta];
  TCanvas *myc1[neta];

  for (unsigned ie(0); ie<neta;++ie){//loop on version
    std::ostringstream lname;
    lname << "myc_" << eta[ie]; 
    myc[ie] = new TCanvas(lname.str().c_str(),lname.str().c_str(),1500,1000);
    lname.str("");
    lname << "myc1_" << eta[ie]; 
    myc1[ie] = new TCanvas(lname.str().c_str(),lname.str().c_str(),1500,1000);
  }
  TGraphErrors *grDummy = new TGraphErrors();
  grDummy->SetMinimum(0);
  grDummy->SetMaximum(1.1);
  TGraphErrors *grDummy1 = new TGraphErrors();
  grDummy1->SetMinimum(isPi ?0 : 0.5);
  grDummy1->SetMaximum(isPi ? 26 : 1000);
  for (unsigned ipt(0);ipt<npt;++ipt){//loop on pt
    grDummy->SetPoint(ipt,pt[ipt],0);
    grDummy1->SetPoint(ipt,pt[ipt],0);
  }
  std::ostringstream label;
  label << ";p_{T} (GeV);proba(n_{E>60fC}#geq" << limit << ")";
  grDummy->SetTitle(label.str().c_str());
  label.str("");
  label << ";p_{T} (GeV);<n_{E>60fC}>";
  grDummy1->SetTitle(label.str().c_str());

  TGraphErrors *gr[nV][neta];
  TGraphErrors *grMean[nV][neta];
  TTree *tree = 0;
  TFile *ftmp = 0;

  TLatex lat;
  char buf[200];

  TLegend *leg = new TLegend(0.68,0.15,0.92,0.4);
  leg->SetFillColor(10);
  for (unsigned ieta(0); ieta<neta;++ieta){//loop on eta
    myc[ieta]->cd();
    grDummy->Draw("AL");
    myc1[ieta]->cd();
    if (!isPi) gPad->SetLogy(1);
    gPad->SetGridy(1);
    grDummy1->Draw("AL");
    grDummy1->GetYaxis()->SetNdivisions(121);
    
    for (unsigned iv(0); iv<nV;++iv){//loop on version
      std::ostringstream lname;
      lname << "grv" << v[iv]<< "_eta" << eta[ieta];
      gr[iv][ieta] = new TGraphErrors();
      gr[iv][ieta]->SetName(lname.str().c_str());
      grMean[iv][ieta] = new TGraphErrors();
      lname << "_mean";
      grMean[iv][ieta]->SetName(lname.str().c_str());
      //set dummy values to have full range for sure...
      for (unsigned ipt(0);ipt<npt;++ipt){//loop on pt
	std::ostringstream fname;
	fname << path << "v" << v[iv] << "_et" << pt[ipt] << "_eta" << eta[ieta];
	if (eta[ieta]==2) fname << ".0";
	fname  << "/NAbove26.root";
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
	if (nentries>0) {
	  std::cout << " mean " << hprob->GetMean() << " " << hprob->GetMeanError() ;
	  grMean[iv][ieta]->SetPoint(grMean[iv][ieta]->GetN(),pt[ipt],hprob->GetMean());
	  grMean[iv][ieta]->SetPointError(grMean[iv][ieta]->GetN()-1,0,hprob->GetMeanError());
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
	gr[iv][ieta]->SetPoint(gr[iv][ieta]->GetN(),pt[ipt],integral);
	gr[iv][ieta]->SetPointError(gr[iv][ieta]->GetN()-1,0,error);
	hprob->Delete();
	ftmp->Close();
      }//loop on pt
      myc[ieta]->cd();
      std::cout << "v " << v[iv] 
		<< " eta " << eta[ieta]
		<< " npoints = " << gr[iv][ieta]->GetN()
		<< " nmean = " << grMean[iv][ieta]->GetN()
		<< std::endl;

      gr[iv][ieta]->SetMarkerColor(iv+1);
      gr[iv][ieta]->SetLineColor(iv+1);
      gr[iv][ieta]->SetMarkerStyle(20+iv);
      gr[iv][ieta]->SetMinimum(0.0);
      gr[iv][ieta]->SetMaximum(1.1);
      gr[iv][ieta]->Draw("PE");

      myc1[ieta]->cd();
      grMean[iv][ieta]->SetMarkerColor(iv+1);
      grMean[iv][ieta]->SetLineColor(iv+1);
      grMean[iv][ieta]->SetMarkerStyle(20+iv);
      grMean[iv][ieta]->SetMinimum(0.01);
      grMean[iv][ieta]->SetMaximum(25);
      grMean[iv][ieta]->Draw("PEL");

    }//loop on versions
    myc[ieta]->cd();
    if (ieta==0) {
      leg->AddEntry(gr[0][ieta],"TP-28-12","P");
      leg->AddEntry(gr[1][ieta],"TP-24-11","P");
      leg->AddEntry(gr[2][ieta],"TP-18-09","P");
    }

    leg->Draw("same");
    sprintf(buf,"#pi^{-}, #eta = %3.2f",eta[ieta]);
    lat.SetTextSize(0.06);
    lat.DrawLatexNDC(0.7,0.5,buf);
    lat.SetTextSize(0.03);
    lat.DrawLatexNDC(0.02,0.02,"HGCAL Standalone Simulation");
    myc[ieta]->Update();
    std::ostringstream lsave;
    lsave << "SummaryNabove60fC_" << particle << "_limit" << limit << "_vspt_eta" << eta[ieta] << ".pdf";
    myc[ieta]->Print(lsave.str().c_str());

    myc1[ieta]->cd();
    leg->Draw("same");
    lat.SetTextSize(0.06);
    lat.DrawLatexNDC(0.15,0.85,buf);
    lat.SetTextSize(0.03);
    lat.DrawLatexNDC(0.02,0.02,"HGCAL Standalone Simulation");
    myc1[ieta]->Update();
    lsave.str("");
    lsave << "SummaryNMeanabove60fC_" << particle << "_vspt_eta" << eta[ieta] << ".pdf";
    myc1[ieta]->Print(lsave.str().c_str());


  }//loop on eta

  return 0;
}//main
