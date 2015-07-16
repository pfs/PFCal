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
#include "TRandom3.h"

double absWeight(const unsigned layer, const bool dedx=true){
  if (layer == 0) return 0.248;
  if (layer == 1) return 1;//95.4/95.4=1
    if (layer == 2) return 0.92;//88.16/95.4=0.92
    if (layer == 3) return 0.537;//51.245/95.4=0.537
    if (layer == 4) return 0.92;
    if (layer == 5) return 0.537;
    if (layer == 6) return 0.92;
    if (layer == 7) return 0.537;
    if (layer == 8) return 0.92;
    if (layer == 9) return 0.537;
    if (layer == 10) return 0.92;
    if (layer == 11) return 0.78;//74.45/95.4=0.78
    if (layer == 12) return 1.071;//102.174/95.4=1.071
    if (layer == 13) return 0.78;
    if (layer == 14) return 1.071;
    if (layer == 15) return 0.78;
    if (layer == 16) return 1.071;
    if (layer == 17) return 0.78;
    if (layer == 18) return 1.071;
    if (layer == 19) return 0.78;
    if (layer == 20) return 1.071;
    if (layer == 21) return 1.1047;//105.39/95.4=1.1047
    if (layer == 22) return 1.378;//131.476/95.4=1.378
    if (layer == 23) return 1.1047;
    if (layer == 24) return 1.378;
    if (layer == 25) return 1.1047;
    if (layer == 26) return 1.378;
    if (layer == 27) return 1.1047;
    if (layer == 28) return 1.378;
    if (layer == 29) return 1.1047;
};
int getNoise(){//main

  TH1F *noiseEE = new TH1F("noiseEE",";E (mips)",100,0,300);
  TH1F *noiseFH = new TH1F("noiseFH",";E (mips)",100,0,300);
  TH1F *noiseBH = new TH1F("noiseBH",";E (mips)",100,0,300);
  TH1F *noise = new TH1F("noise",";E (mips)",100,0,300);
  TH1F *aboveEE = new TH1F("aboveEE",";E (mips)",100,0,1000);
  TH1F *aboveFH = new TH1F("aboveFH",";E (mips)",100,0,1000);
  TH1F *aboveBH = new TH1F("aboveBH",";E (mips)",100,0,1000);

  const unsigned nEvts = 1500;
  const unsigned nEE = 100*100;//1*1 m^2
  const unsigned nFH = 100*100*9;
  const unsigned nBH = 50*50*12;

  const double sigSi = 0.14;
  const double sigScint = 0.2;

  TRandom3 lRand;
  lRand.SetSeed(0);

  for (unsigned ie(0); ie<nEvts;++ie){//loop on events
    double etot_ee = 0;
    unsigned nAbove = 0;
    for (unsigned iL(0); iL<30;++iL){
      for (unsigned ee(0); ee<nEE;++ee){
	double e = lRand.Gaus(0,sigSi);
	if (e>0.5) {
	  etot_ee += e*absWeight(iL);
	  nAbove++;
	}
      }
    }
    noiseEE->Fill(etot_ee);
    aboveEE->Fill(nAbove);
    double etot_fh = 0;
    nAbove = 0;
    for (unsigned fh(0); fh<nFH;++fh){
      double e = lRand.Gaus(0,sigSi);
      if (e>0.5) {
	etot_fh += e;
	nAbove++;
      }
    }
    noiseFH->Fill(etot_fh);
    aboveFH->Fill(nAbove);
    double etot_bh = 0;
    nAbove = 0;
    for (unsigned bh(0); bh<nBH;++bh){
      double e = lRand.Gaus(0,sigScint);
      if (e>0.5) {
	etot_bh += e;
	nAbove++;
      }
    }
    noiseBH->Fill(etot_bh);
    aboveBH->Fill(nAbove);
    noise->Fill(etot_bh+etot_fh+etot_ee);

  }//loop on events

  noise->Draw();
  noiseEE->SetLineColor(4);
  noiseEE->Draw("same");
  noiseFH->SetLineColor(2);
  noiseFH->Draw("same");
  noiseBH->SetLineColor(3);
  noiseBH->Draw("same");

  aboveEE->SetLineStyle(3);
  aboveEE->SetLineColor(4);
  aboveEE->Draw("same");
  aboveFH->SetLineStyle(3);
  aboveFH->SetLineColor(2);
  aboveFH->Draw("same");
  aboveBH->SetLineStyle(3);
  aboveBH->SetLineColor(3);
  aboveBH->Draw("same");


  std::cout << "Det & Mean & RMS \\\\" << std::endl
	    << "EE & " << noiseEE->GetMean() << "&" << noiseEE->GetRMS() << "\\\\" << aboveEE->GetMean() << std::endl
	    << "FH & " << noiseFH->GetMean() << "&" << noiseFH->GetRMS() << "\\\\" << aboveFH->GetMean() << std::endl
	    << "BH & " << noiseBH->GetMean() << "&" << noiseBH->GetRMS() << "\\\\" << aboveBH->GetMean() << std::endl
	    << "All & " << noise->GetMean() << "&" << noise->GetRMS() << "\\\\" << std::endl;


  return 0;
}//main
