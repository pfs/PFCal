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
#include "TProfile.h"
#include "effSigmaMacro.C"

double absWeight(const unsigned layer, bool dedx=true){
  if (!dedx){
    if (layer == 0) return 0.08696;
    if (layer == 1) return 1;//0.92
    if (layer == 2) return 0.646989;//88.16/95.4=0.92
    if (layer == 3) return 0.617619;//51.245/95.4=0.537
    if (layer == 4) return 0.646989;
    if (layer == 5) return 0.617619;
    if (layer == 6) return 0.646989;
    if (layer == 7) return 0.617619;
    if (layer == 8) return 0.646989;
    if (layer == 9) return 0.617619;
    if (layer == 10) return 0.646989;
    if (layer == 11) return 0.942829;//74.45/95.4=0.78
    if (layer == 12) return 0.859702;//102.174/95.4=1.071
    if (layer == 13) return 0.942829;
    if (layer == 14) return 0.859702;
    if (layer == 15) return 0.942829;
    if (layer == 16) return 0.859702;
    if (layer == 17) return 0.942829;
    if (layer == 18) return 0.859702;
    if (layer == 19) return 0.942829;
    if (layer == 20) return 0.859702;
    if (layer == 21) return 1.37644;//105.39/95.4=1.1047
    if (layer == 22) return 1.30447;//131.476/95.4=1.378
    if (layer == 23) return 1.37644;
    if (layer == 24) return 1.30447;
    if (layer == 25) return 1.37644;
    if (layer == 26) return 1.30447;
    if (layer == 27) return 1.37644;
    if (layer == 28) return 1.30447;
    if (layer == 29) return 1.37644;//1.79662;//
  } else {
    if (layer == 0) return 0.248;//0.06588;
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
  }
  return 1;
};


void removeOneLayer(TCanvas **myc,
		    TGraphErrors* & gr,
		    TTree *t,
		    const unsigned nLayers,
		    const std::vector<unsigned>& lToRemove,
		    unsigned & min){

  TLatex lat;
  char buf[100];
  unsigned nP = 31;
  double layer[nP];
  double reso[nP];
  double layerErr[nP];
  double resoErr[nP];

  bool doG1 = true;
  bool doG2 = true;

  std::string cut1 = "converted1==0 && eTrue1/cosh(etaTrue1)>40 && TMath::Abs(etaTrue1)>1.6 && TMath::Abs(etaTrue1)<2.8 && (TMath::Nint(TMath::Abs(phiTrue1)/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(phiTrue1)/(2*TMath::Pi())*360.)%20>11) && invalidDetid1==0 && invalidNeighbour1==0";
  std::string cut2 = "converted2==0 && eTrue2/cosh(etaTrue2)>40 && TMath::Abs(etaTrue2)>1.6 && TMath::Abs(etaTrue2)<2.8 && (TMath::Nint(TMath::Abs(phiTrue2)/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(phiTrue2)/(2*TMath::Pi())*360.)%20>11) && invalidDetid2==0 && invalidNeighbour2==0";

  std::ostringstream label;
  label.str("");
  label << ";e_{3#times3};photons";

  myc[7]->cd();
  std::ostringstream histname;
  histname << "calib" << lToRemove.size() ;
  TH2F *calib = new TH2F(histname.str().c_str(),";e_{gen};e_{reco}",100,10,500,100,1000,40000);

  std::ostringstream label1tot;
  std::ostringstream label2tot;
  label1tot.str("");
  label1tot << "(";
  label2tot.str("");
  label2tot << "(";
  bool firstLayer = true;
  for (unsigned jL(0);jL<nLayers;++jL){
    if (!firstLayer) {
      label1tot << "+";
      label2tot << "+";
    }
    label1tot << "eSR3RecoLayer"<<jL<<"_1*absWeight(" << jL << ")";
    label2tot << "eSR3RecoLayer"<<jL<<"_2*absWeight(" << jL << ")";
    firstLayer=false;
  }//loop on layers
  label1tot << ")";///eTrue1";
  label2tot << ")";///eTrue2";

  t->Draw((label1tot.str()+":eTrue1>>"+histname.str()).c_str(),cut1.c_str());
 
  //t->Draw(("eSR3Reco1:eTrue1>>"+histname.str()).c_str(),cut1.c_str());
  //TH1F *hist = (TH1F*)(gPad->GetPrimitive("htemp"))->Clone(histname.str().c_str());
  //hist->SetTitle(label.str().c_str());
  t->Draw((label2tot.str()+":eTrue2>>+"+histname.str()).c_str(),cut2.c_str());
    //t->Draw(("eSR3Reco2:eTrue2>>+"+histname.str()).c_str(),cut2.c_str());
  calib->Draw("colz");
  calib->ProfileX();
  gPad->Update();
  histname << "_pfx";
  TProfile *h2_pfx = (TProfile*)gDirectory->Get(histname.str().c_str());
  h2_pfx->SetMarkerStyle(22);
  h2_pfx->SetMarkerColor(1);
  h2_pfx->Draw("PEsame");
  h2_pfx->Fit("pol1","RIQ","same",100,400);

  TF1 *fit = h2_pfx->GetFunction("pol1");
  double slopeTot = fit->GetParameter(1);
  double offsetTot = fit->GetParameter(0);

  myc[5]->cd();
  histname.str("");
  histname << "totalE" << lToRemove.size() ;
  TH1F *totalE = new TH1F(histname.str().c_str(),";e_{3#times3};photons",150,0.8,1.2);
 std:ostringstream var;

  var << "(" << label1tot.str() << "-" << offsetTot << ")/(" << slopeTot << "*eTrue1)>>" << histname.str();

  //var << "(eSR3Reco1-" << offsetTot << ")/(" << slopeTot << "*eTrue1)>>" << histname.str();
  if (doG1) t->Draw(var.str().c_str(),cut1.c_str());
  //TH1F *hist = (TH1F*)(gPad->GetPrimitive("htemp"))->Clone(histname.str().c_str());
  //hist->SetTitle(label.str().c_str());
  var.str("");
  var << "(" << label2tot.str() << "-" << offsetTot << ")/(" << slopeTot << "*eTrue2)>>" << histname.str();
  //var << "(eSR3Reco2-" << offsetTot << ")/(" << slopeTot << "*eTrue2)>>+" << histname.str();
  if (doG2) t->Draw(var.str().c_str(),cut2.c_str());
  totalE->Draw("PE");

  totalE->Fit("gaus","Q");//,"",
  //totalE->GetMean()-1.5*totalE->GetRMS(),
  //	      totalE->GetMean()+1.5*totalE->GetRMS());
  fit = totalE->GetFunction("gaus");
  //std::cout << " -- Starting point: " << fit->GetParameter(2)/fit->GetParameter(1)*100
  //<< " \\pm " <<  fit->GetParError(2)/sqrt(2*totalE->GetEntries())*100
  //<< std::endl;
  layer[0] = -1;
  layerErr[0] = 0;
  reso[0] = fit->GetParameter(2)/fit->GetParameter(1)*100;
  resoErr[0] = fit->GetParameter(2)/sqrt(2*totalE->GetEntries())*100;

  std::ostringstream label1[nLayers];
  std::ostringstream label2[nLayers];
  for (unsigned iL(0);iL<nLayers;++iL){
    label1[iL].str("");
    label1[iL] << "(";
    label2[iL].str("");
    label2[iL] << "(";
    firstLayer = true;
    std::ostringstream absSum;
    for (unsigned jL(0);jL<nLayers;++jL){
      if (jL==iL) {
	absSum << "absWeight(" << jL << ")+";
	continue;
      }
      bool skipLayer = false;
      for (unsigned r(0); r<lToRemove.size(); r++){
	if (jL==lToRemove[r]) {
	  skipLayer = true;
	  absSum << "absWeight(" << jL << ")+";
	  break;
	}
      }
      if (skipLayer) continue;
      if (!firstLayer) {
	label1[iL] << "+";
	label2[iL] << "+";
      }
      label1[iL] << "eSR3RecoLayer"<<jL<<"_1*(" << absSum.str() << "absWeight(" << jL << "))";
      label2[iL] << "eSR3RecoLayer"<<jL<<"_2*(" << absSum.str() << "absWeight(" << jL << "))";
      absSum.str("");
      firstLayer=false;
    }//loop on layers
    label1[iL] << ")";///eTrue1";
    label2[iL] << ")";///eTrue2";
    //std::cout << " Layer " << iL << ": "
    //<< label1[iL].str() << std::endl;
  }//loop on layers
  TH1F *hist[nLayers];
  TH2F *calibr[nLayers];
  min = 1;
  double minval = 1000;
  for (unsigned iL(0);iL<nLayers;++iL){
    myc[iL/6]->cd(iL%6+1);
    histname.str("");
    histname << "calib" << lToRemove.size() << "_" << iL ;
    calibr[iL] = new TH2F(histname.str().c_str(),";e_{gen};e_{reco}",100,10,500,100,1000,40000);
    t->Draw((label1[iL].str()+":eTrue1>>"+histname.str()).c_str(),cut1.c_str());
    //TH1F *hist = (TH1F*)(gPad->GetPrimitive("htemp"))->Clone(histname.str().c_str());
    //hist->SetTitle(label.str().c_str());
    t->Draw((label2[iL].str()+":eTrue2>>+"+histname.str()).c_str(),cut2.c_str());
    calibr[iL]->Draw("colz");
    calibr[iL]->ProfileX();
    gPad->Update();
    histname << "_pfx";
    h2_pfx = (TProfile*)gDirectory->Get(histname.str().c_str());
    h2_pfx->SetMarkerStyle(22);
    h2_pfx->SetMarkerColor(1);
    h2_pfx->Draw("PEsame");
    h2_pfx->Fit("pol1","RIQ","same",100,500);
    gPad->Update();
    fit = h2_pfx->GetFunction("pol1");
    double slope = fit->GetParameter(1);
    double offset = fit->GetParameter(0);
    
    
    myc[iL/6]->cd(iL%6+1);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    histname.str("");
    histname << "hist" << lToRemove.size() << "_" << iL;
    hist[iL] = new TH1F(histname.str().c_str(),";e_{3#times3}/e_{true};photons",150,0.8,1.2);
    var.str("");
    var << "(" << label1[iL].str() << "-" << offset << ")/(" << slope << "*eTrue1)>>" << histname.str();
    //var << "(eSR3Reco1-" << label1[iL].str() << ")/eSR3Reco1>>" << histname.str();
    if (doG1) t->Draw(var.str().c_str(),cut1.c_str());
    //hist = (TH1F*)(gPad->GetPrimitive("htemp"))->Clone(histname.str().c_str());
    //std::cout << " photon 1 " << hist[iL]->GetEntries() << std::endl;
    var.str("");
    var << "(" << label2[iL].str() << "-" << offset << ")/(" << slope << "*eTrue2)>>+" << histname.str();
    //var << "(eSR3Reco2-" << label2[iL].str() << ")/eSR3Reco2>>" << histname.str();
    if (doG2) t->Draw(var.str().c_str(),cut2.c_str());
    hist[iL]->SetTitle(label.str().c_str());
    hist[iL]->Draw();
    //std::cout << " photons 1+2 " << hist[iL]->GetEntries() << std::endl;
    hist[iL]->Fit("gaus","Q");//,"",
    //hist[iL]->GetMean()-1.5*hist[iL]->GetRMS(),
    //		  hist[iL]->GetMean()+1.5*hist[iL]->GetRMS());
    fit = hist[iL]->GetFunction("gaus");
    //std::cout << " -- Layer " << iL << " " << fit->GetParameter(2)/fit->GetParameter(1)*100
    // << " \\pm " <<  fit->GetParError(2)/sqrt(2*hist[iL]->GetEntries())*100
    //<< std::endl;
    layer[iL+1] = iL;
    layerErr[iL+1] = 0;
    reso[iL+1] = fit->GetParameter(2)*100;///fit->GetParameter(1)*100;
    resoErr[iL+1] = fit->GetParError(2)*100;//ameter(2)/sqrt(2*hist[iL]->GetEntries())*100;

    if (iL>0 && iL<nLayers-1 && reso[iL+1]<minval){
      bool ignore = false;
      for (unsigned r(0); r<lToRemove.size(); r++){
	if (abs(static_cast<int>(iL)-static_cast<int>(lToRemove[r]))<2) ignore = true;
      }
      if (!ignore){
	min = iL;
	minval = reso[iL+1];
      }
    }

    sprintf(buf,"Layer %d",iL);
    lat.DrawLatex(hist[iL]->GetXaxis()->GetBinCenter(hist[iL]->GetNbinsX()/2),0.9*hist[iL]->GetMaximum(),buf);
    myc[iL/6]->Update();
  }

  gr = new TGraphErrors(nP,layer,reso,layerErr,resoErr);
  gr->SetTitle(";Layer removed;#sigma(#DeltaE/E)");
  gr->SetMinimum(1.5);
  gr->SetMaximum(4);
  myc[6]->cd();
  //gr->Draw("AP");

};


int plotDescopingChoice(){
  TFile *fin = TFile::Open("PCA_topoFix_Hgg_0pu.root");
  if (!fin) return 1;
  fin->cd("hgg");
  const unsigned nLayers = 30;

  TTree *t = (TTree*)gDirectory->Get("tree");
  
  if (!t) {
    std::cout << " Error, tree not found!" << std::endl;
    return 1;
  }
  
  gStyle->SetOptStat("eMRuo");
  gStyle->SetOptFit(1111);


  const unsigned nC = 8;
  TCanvas *myc[nC];
  for (unsigned ic(0); ic<nC; ++ic){
    std::ostringstream label;
    label << "myc" << ic ;
    myc[ic] = new TCanvas(label.str().c_str(),
			  label.str().c_str(),
			  1500,1000);
    if (ic<5) myc[ic]->Divide(3,2);
  }

  TGraphErrors *gr1 = 0;
  std::vector<unsigned> lToRemove;
  unsigned inimin;
  removeOneLayer(myc,gr1,t,nLayers,lToRemove,inimin);
  gr1->SetName("grRef");
  TLatex lat;
  char buf[100];

  myc[6]->cd();
  gPad->SetGridy(1);
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerColor(1);
  gr1->SetLineColor(1);
  gr1->Draw("AP");
  myc[6]->Update();


  const unsigned nRemove = 12;
  //unsigned list[nRemove] = {22,24,26,28,3,10,1,16,14,5,7,18};//,5,20,7,9,11,18};
  TGraphErrors *gr[nRemove];

  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->SetFillStyle(0);
  leg->AddEntry(gr1,"All layers","P");


  lToRemove.push_back(inimin);

  std::cout << " -- -1 Layer min: " << inimin << std::endl;

  std::ostringstream label;
  unsigned r=0;
  for (; r<nRemove; r++){
    label << "-L" << lToRemove[r];
    //}
    //r--;
    gr[r] = 0;
    unsigned newmin = 1;
    removeOneLayer(myc,gr[r],t,nLayers,lToRemove,newmin);
    std::ostringstream name;
    name << "gr" << r;
    gr[r]->SetName(name.str().c_str());
    
    myc[6]->cd();
    gr[r]->SetMarkerStyle(21+r);
    if (r%7!=3){
      gr[r]->SetMarkerColor(2+r%7);
      gr[r]->SetLineColor(2+r%7);
    }
    else {
      gr[r]->SetMarkerColor(9);
      gr[r]->SetLineColor(9);
    }
    gr[r]->Draw("Psame");
    myc[6]->Update();
    //label << "-L" << list[r];
    leg->AddEntry(gr[r],label.str().c_str(),"P");
    lToRemove.push_back(newmin);
    std::cout << " -- -" << r+2 << " Layer min: " << newmin << std::endl;
  }

  leg->Draw("same");

  myc[6]->Update();
  myc[6]->Print("Descop/Scenarios_reso_auto_all_dedx.pdf");

  /*myc[7]->cd();
  TGraphErrors *gr2[nLayers];
  //unsigned iL(0);
  for (unsigned iL(0);iL<nLayers;++iL){
  //for (unsigned iL(0);iL<2;++iL){
    std::cout << " -- Removing layer " << iL << std::endl;
    gr2[iL] = 0;
    lToRemove.clear();
    lToRemove.push_back(iL);
    removeOneLayer(myc,gr2[iL],t,nLayers,lToRemove);
    std::ostringstream name;
    name << "gr" << iL;
    gr2[iL]->SetName(name.str().c_str());
    gr2[iL]->SetMarkerStyle(20+iL/9);
    gr2[iL]->SetMarkerColor(1+iL%9);
    gr2[iL]->SetLineColor(1+iL%9);
    gr2[iL]->SetMinimum(1);
    gr2[iL]->SetMaximum(4);
    gr2[iL]->Draw(iL==0?"AP":"Psame");
    //return 1;
    }
  myc[7]->Update();
  myc[7]->Print("Descop/TwoLayers.pdf");
  */
  //for (unsigned ic(0); ic<nC; ++ic){
  //std::ostringstream label;
    //label << "PLOTS/LatProf_Hgg_0pu_" << numSR << "over" << denSR << "_l"<<6*ic << "-" << 6*(ic+1)-1 << ".pdf";
    //myc[ic]->Print(label.str().c_str());
  //}

  return 0;
}
