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
    if (layer == 0) return 0.12;//0.06588;
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

double calibratedE(const double Etot, const double eta){
  //calibration for signal region 2: 3*3 cm^2
  double pars[3] = {69.5,4.5,-0.8};
  //double pars[3] = {75.5,0,0};
  double paro[3] = {-34.4,0,0};
  //double pars[3] = {77,3.4,-0.50};
  //double paro[3] = {-11.6,-7.7,-8.8};
  //double paro[3] = {-5.3,-12.8,-6.9};  
  double offset = paro[0] + paro[1]*fabs(eta) + paro[2]*eta*eta;
  double slope = pars[0] + pars[1]*fabs(eta) + pars[2]*eta*eta;
  return (Etot-offset)/slope;
};

int plotDescopingTree(){
  bool doCMSSW = true;
  //bool doDEDX = true;
  std::string savestr = doCMSSW? "Hgg" : "eta21_3x3";
  //if (!doDEDX) savestr += "_x0";
  //else savestr += "_dedx";
  //savestr += "_inverseList";
  savestr += "_dedxchoice";

  for (unsigned iL(0); iL<30; ++iL){
    std::cout << " " << iL;
    if (iL<10) std::cout << " ";
    std::cout << " & " << absWeight(iL,false) << " & " << absWeight(iL,true) << " \\\\" << std::endl;
  }

  double eta = 2.1;
  TFile *fin;
  if (doCMSSW) fin = TFile::Open("PCA_topoFix_Hgg_0pu.root");
  else fin = TFile::Open("eta21_pu0.root");

  const unsigned nE = 17;
  double etrue[nE] = {100,10,125,150,175,200,20,30,3,40,50,5,60,70,7,80,90};

  if (!fin) return 1;
  if (doCMSSW) fin->cd("hgg");
  else fin->cd("Energies");

  const unsigned nLayers = 30;

  TTree *t;
  if (doCMSSW) t = (TTree*)gDirectory->Get("tree");
  else t = (TTree*)gDirectory->Get("Ereso");
  if (!t) {
    std::cout << " Error, tree not found!" << std::endl;
    return 1;
  }
  
  gStyle->SetOptStat("eMRuo");
  gStyle->SetOptFit(1111);
  gStyle->SetStatW(0.3);


  std::ostringstream label;

  const unsigned nC = 4;
  TCanvas *myc[nC];
  for (unsigned ic(0); ic<nC; ++ic){
    label.str("");
    label << "myc" << ic ;
    myc[ic] = new TCanvas(label.str().c_str(),
			  label.str().c_str(),
			  1800,1000);
    if (ic<3) myc[ic]->Divide(5,3);
    //else myc[ic]->Divide(2,1);
  }


  ////////////////////////
  //do calibration
  ////////////////////////

  label.str("");
  label << ";e_{3#times3};photons";

  const unsigned nRemove = 12;
  //x0
  //unsigned list[nRemove] = {22,7,5,28,1,3,26,9,18,16,20,24};
  //dedx
  unsigned list[nRemove] = {25,27,15,1,10,3,18,5,12,7,23,20};
  unsigned color[12]={2,3,4,6,7,8,9,30,38,40,42,50};
  TH1F *histx0[nRemove+1];
  TH1F *histdedx[nRemove+1];
  //TH1F *histE[nRemove+1][nE];
  TH1F *histcalibx0[nRemove+1];
  TH2F *calibrx0[nRemove+1];
  TH1F *histcalibdedx[nRemove+1];
  TH2F *calibrdedx[nRemove+1];
  std::ostringstream histname;
  for (unsigned ir(0); ir<nRemove+1; ir++){
    histname.str("");
    histname << "calibx0" << "_" << ir ;
    calibrx0[ir] = new TH2F(histname.str().c_str(),";e_{gen};e_{reco}",100,10,500,100,1000,40000);
    histname.str("");
    histname << "histx0" << "_" << ir;
    histx0[ir] = new TH1F(histname.str().c_str(),";e_{rec}/eTrue;photons",300,50,150);
    histname.str("");
    histname << "histCalibx0" << "_" << ir;
    histcalibx0[ir] = new TH1F(histname.str().c_str(),";e_{rec}/eTrue;photons",300,0.8,1.6);
    histname.str("");
    histname << "calibdedx" << "_" << ir ;
    calibrdedx[ir] = new TH2F(histname.str().c_str(),";e_{gen};e_{reco}",100,10,500,100,1000,40000);
    histname.str("");
    histname << "histdedx" << "_" << ir;
    histdedx[ir] = new TH1F(histname.str().c_str(),";e_{rec}/eTrue;photons",300,50,150);
    histname.str("");
    histname << "histCalibdedx" << "_" << ir;
    histcalibdedx[ir] = new TH1F(histname.str().c_str(),";e_{rec}/eTrue;photons",300,0.8,1.6);
    //for (unsigned ie(0); ie<nE;++ie){
    //histname.str("");
    //histname << "hist" << "_" << ir << "_e" << etrue[ie];
    //histE[ir][ie] = new TH1F(histname.str().c_str(),";e_{rec};photons",1000,100,40000);
    //}
  }
  
  const unsigned nPhotons = doCMSSW ? 2 : 1;

  float eRecoLayers[nPhotons][nLayers];
  unsigned eventIndex;
  float etaTrue[nPhotons],phiTrue[nPhotons],eTrue[nPhotons];
  int converted[nPhotons];
  char invalidDetid[nPhotons],invalidNeighbour[nPhotons];
  if (doCMSSW){
    for (unsigned ip(0); ip<nPhotons;++ip){
      std::string suffix = ip==0?"1" : "2";
      t->SetBranchAddress(("etaTrue"+suffix).c_str(),&etaTrue[ip]);
      t->SetBranchAddress(("eTrue"+suffix).c_str()  ,&eTrue[ip]);
      t->SetBranchAddress(("phiTrue"+suffix).c_str()  ,&phiTrue[ip]);
      t->SetBranchAddress(("converted"+suffix).c_str(),&converted[ip]);
      t->SetBranchAddress(("invalidDetid"+suffix).c_str(),&invalidDetid[ip]);
      t->SetBranchAddress(("invalidNeighbour"+suffix).c_str(),&invalidNeighbour[ip]);
      
      for (unsigned iL(0); iL<nLayers;++iL){
	label.str("");
	label << "eSR3RecoLayer"<<iL<<"_"<<suffix;
	t->SetBranchAddress(label.str().c_str(),&eRecoLayers[ip][iL]);
      }
    }
  }
  else {
    t->SetBranchAddress("eventIndex",&eventIndex);
    for (unsigned iL(0); iL<nLayers;++iL){
      label.str("");
      label << "energy_"<<iL<<"_SR2";
      t->SetBranchAddress(label.str().c_str(),&eRecoLayers[0][iL]);
    }
  }

  //fill histograms
  const unsigned nEvts = t->GetEntries();
  unsigned eIdx = 0;
  for (unsigned ievt(0); ievt<nEvts;++ievt){//loop on events
    t->GetEntry(ievt);

    bool pass[nPhotons];
    if (doCMSSW){
      pass[0] = converted[0]==0 && eTrue[0]/cosh(etaTrue[0])>40 && TMath::Abs(etaTrue[0])>1.6 && TMath::Abs(etaTrue[0])<2.8 && (TMath::Nint(TMath::Abs(phiTrue[0])/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(phiTrue[0])/(2*TMath::Pi())*360.)%20>11) && invalidDetid[0]==0 && invalidNeighbour[0]==0;
      pass[1] = converted[1]==0 && eTrue[1]/cosh(etaTrue[1])>40 && TMath::Abs(etaTrue[1])>1.6 && TMath::Abs(etaTrue[1])<2.8 && (TMath::Nint(TMath::Abs(phiTrue[1])/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(phiTrue[1])/(2*TMath::Pi())*360.)%20>11) && invalidDetid[1]==0 && invalidNeighbour[1]==0;
    }
    else {
      //std::cout << ievt << " " << eventIndex << " " << etrue[eIdx] << std::endl;
      if (eIdx<nE) eTrue[0] = etrue[eIdx]*cosh(eta);
      else return 1;
      pass[0] = true;
      if (ievt>0 && eventIndex==0) eIdx++;
    }

    for (unsigned ip(0); ip<nPhotons;++ip){//loop on photons
      if (!pass[ip]) continue;
      for (unsigned ir(0); ir<nRemove+1; ir++){
	double etotx0 = 0;
	double etotdedx = 0;
	double absSumx0 = 0;
	double absSumdedx = 0;
	for (unsigned jL(0);jL<nLayers;++jL){
	  if (ir>0){
	    bool skipLayer = false;
	    for (unsigned r(0); r<ir; r++){
	      if (jL==list[r]) {
		skipLayer = true;
		absSumx0+=absWeight(jL,false);
		absSumdedx+=absWeight(jL,true);
		break;
	      }
	    }
	    if (skipLayer) continue;
	  }
	  absSumx0+=absWeight(jL,false);
	  absSumdedx+=absWeight(jL,true);
	  etotx0 += eRecoLayers[ip][jL]*absSumx0;
	  etotdedx += eRecoLayers[ip][jL]*absSumdedx;
	  absSumx0=0;
	  absSumdedx=0;
	}//loop on layers
	histx0[ir]->Fill(etotx0/eTrue[ip]);
	calibrx0[ir]->Fill(eTrue[ip],etotx0);
 	histdedx[ir]->Fill(etotdedx/eTrue[ip]);
	calibrdedx[ir]->Fill(eTrue[ip],etotdedx);
     }//loop on remove
    }//loop on photons

  }//loop on events

  double slope[2][nRemove+1];
  double offset[2][nRemove+1];
  double slopeerr[2][nRemove+1];
  double offseterr[2][nRemove+1];
  double reso[2][nRemove+1];
  double resoerr[2][nRemove+1];
  double resoCalib[2][nRemove+1];
  double resoCaliberr[2][nRemove+1];
  double nDropped[nRemove+1];
  double nDroppederr[nRemove+1];

  for (unsigned ir(0); ir<nRemove+1; ir++){
    myc[0]->cd(ir+1);
    calibrx0[ir]->Draw("colz");
    calibrx0[ir]->ProfileX();
    gPad->Update();
    histname.str("");
    histname << "calibx0" << "_" << ir << "_pfx";
    TProfile *h2_pfx = (TProfile*)gDirectory->Get(histname.str().c_str());
    h2_pfx->SetMarkerStyle(22);
    h2_pfx->SetMarkerColor(1);
    h2_pfx->Draw("PEsame");
    if (!doCMSSW) h2_pfx->Fit("pol1","RIQ","same",10,400);
    else h2_pfx->Fit("pol1","RIQ","same",100,500);
    gPad->Update();
    TF1 *fit = h2_pfx->GetFunction("pol1");
    slope[0][ir] = fit->GetParameter(1);
    slopeerr[0][ir] = fit->GetParError(1);
    offset[0][ir] = fit->GetParameter(0);
    offseterr[0][ir] = fit->GetParError(0);
    nDropped[ir] = ir;
    nDroppederr[ir] = 0;

    calibrdedx[ir]->Draw("colz");
    calibrdedx[ir]->ProfileX();
    gPad->Update();
    histname.str("");
    histname << "calibdedx" << "_" << ir << "_pfx";
    h2_pfx = (TProfile*)gDirectory->Get(histname.str().c_str());
    h2_pfx->SetMarkerStyle(22);
    h2_pfx->SetMarkerColor(1);
    h2_pfx->Draw("PEsame");
    if (!doCMSSW) h2_pfx->Fit("pol1","RIQ","same",10,400);
    else h2_pfx->Fit("pol1","RIQ","same",100,500);
    gPad->Update();
    fit = h2_pfx->GetFunction("pol1");
    slope[1][ir] = fit->GetParameter(1);
    slopeerr[1][ir] = fit->GetParError(1);
    offset[1][ir] = fit->GetParameter(0);
    offseterr[1][ir] = fit->GetParError(0);
    nDropped[ir] = ir;
    nDroppederr[ir] = 0;

    myc[1]->cd(ir+1);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    histx0[ir]->Draw();
    histx0[ir]->Fit("gaus","Q","",
		  histx0[ir]->GetMean()-1.5*histx0[ir]->GetRMS(),
		  histx0[ir]->GetMean()+1.5*histx0[ir]->GetRMS());
    fit = histx0[ir]->GetFunction("gaus");
    reso[0][ir] = fit->GetParameter(2)/fit->GetParameter(1);
    resoerr[0][ir] = reso[0][ir]/sqrt(2*histx0[ir]->GetEntries());

    histdedx[ir]->SetLineColor(2);
    histdedx[ir]->Draw("same");
    histdedx[ir]->Fit("gaus","Q","",
		  histdedx[ir]->GetMean()-1.5*histdedx[ir]->GetRMS(),
		  histdedx[ir]->GetMean()+1.5*histdedx[ir]->GetRMS());
    fit = histdedx[ir]->GetFunction("gaus");
    reso[1][ir] = fit->GetParameter(2)/fit->GetParameter(1);
    resoerr[1][ir] = reso[1][ir]/sqrt(2*histdedx[ir]->GetEntries());

  }

  //if (doCMSSW){
  eIdx = 0;
    for (unsigned ievt(0); ievt<nEvts;++ievt){//loop on events
      t->GetEntry(ievt);
      
      bool pass[nPhotons];
      if (doCMSSW){
	pass[0] = converted[0]==0 && eTrue[0]/cosh(etaTrue[0])>40 && TMath::Abs(etaTrue[0])>1.6 && TMath::Abs(etaTrue[0])<2.8 && (TMath::Nint(TMath::Abs(phiTrue[0])/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(phiTrue[0])/(2*TMath::Pi())*360.)%20>11) && invalidDetid[0]==0 && invalidNeighbour[0]==0;
	pass[1] = converted[1]==0 && eTrue[1]/cosh(etaTrue[1])>40 && TMath::Abs(etaTrue[1])>1.6 && TMath::Abs(etaTrue[1])<2.8 && (TMath::Nint(TMath::Abs(phiTrue[1])/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(phiTrue[1])/(2*TMath::Pi())*360.)%20>11) && invalidDetid[1]==0 && invalidNeighbour[1]==0;
      } else {
	if (eIdx<nE) eTrue[0] = etrue[eIdx]*cosh(eta);
	else return 1;
	pass[0] = true;
	if (ievt>0 && eventIndex==0) eIdx++;
      }
      
      for (unsigned ip(0); ip<nPhotons;++ip){//loop on photons
	if (!pass[ip]) continue;
	for (unsigned ir(0); ir<nRemove+1; ir++){
	  double etotx0 = 0;
	  double etotdedx = 0;
	  double absSumx0 = 0;
	  double absSumdedx = 0;
	  for (unsigned jL(0);jL<nLayers;++jL){
	    if (ir>0){
	      bool skipLayer = false;
	      for (unsigned r(0); r<ir; r++){
		if (jL==list[r]) {
		  skipLayer = true;
		  absSumx0+=absWeight(jL,false);
		  absSumdedx+=absWeight(jL,true);
		  break;
		}
	      }
	      if (skipLayer) continue;
	    }
	    absSumx0+=absWeight(jL,false);
	    absSumdedx+=absWeight(jL,true);
	    etotx0 += eRecoLayers[ip][jL]*absSumx0;
	    etotdedx += eRecoLayers[ip][jL]*absSumdedx;
	    absSumx0=0;
	    absSumdedx=0;
	  }//loop on layers
	  if (doCMSSW) histcalibx0[ir]->Fill((etotx0-offset[0][ir])/(slope[0][ir]*eTrue[ip]));
	  else histcalibx0[ir]->Fill(calibratedE(etotx0,eta)/eTrue[ip]);
	  if (doCMSSW) histcalibdedx[ir]->Fill((etotdedx-offset[1][ir])/(slope[1][ir]*eTrue[ip]));
	  //else histcalib[ir]->Fill(calibratedE(etot,eta)/eTrue[ip]);
	}//loop on remove
      }//loop on photons
      
    }//loop on events
    
    
    for (unsigned ir(0); ir<nRemove+1; ir++){
      myc[2]->cd(ir+1);
      gPad->SetGridx(1);
      gPad->SetGridy(1);
      histcalibx0[ir]->Draw();
      histcalibx0[ir]->Fit("gaus","Q");
      TF1 *fit = histcalibx0[ir]->GetFunction("gaus");
      resoCalib[0][ir] = fit->GetParameter(2)/fit->GetParameter(1);
      resoCaliberr[0][ir] = resoCalib[0][ir]/sqrt(2*histcalibx0[ir]->GetEntries());
      histcalibdedx[ir]->Draw();
      histcalibdedx[ir]->Fit("gaus","Q");
      fit = histcalibdedx[ir]->GetFunction("gaus");
      resoCalib[1][ir] = fit->GetParameter(2)/fit->GetParameter(1);
      resoCaliberr[1][ir] = resoCalib[1][ir]/sqrt(2*histcalibdedx[ir]->GetEntries());
    }
    // }


  TGraphErrors *grSlopex0 = new TGraphErrors(nRemove+1,nDropped,slope[0],nDroppederr,slopeerr[0]);
  TGraphErrors *grOffsetx0 = new TGraphErrors(nRemove+1,nDropped,offset[0],nDroppederr,offseterr[0]);
  TGraphErrors *grResox0 = new TGraphErrors(nRemove+1,nDropped,reso[0],nDroppederr,resoerr[0]);
  TGraphErrors *grResoCalibx0 = new TGraphErrors(nRemove+1,nDropped,resoCalib[0],nDroppederr,resoCaliberr[0]);
  TGraphErrors *grSlopededx = new TGraphErrors(nRemove+1,nDropped,slope[1],nDroppederr,slopeerr[1]);
  TGraphErrors *grOffsetdedx = new TGraphErrors(nRemove+1,nDropped,offset[1],nDroppederr,offseterr[1]);
  TGraphErrors *grResodedx = new TGraphErrors(nRemove+1,nDropped,reso[1],nDroppederr,resoerr[1]);
  TGraphErrors *grResoCalibdedx = new TGraphErrors(nRemove+1,nDropped,resoCalib[1],nDroppederr,resoCaliberr[1]);
  myc[0]->cd(14);
  grSlopex0->SetTitle(";# layers removed;slope (mips/GeV)");
  grSlopex0->SetMarkerStyle(20);
  grSlopex0->Draw("AP");
  grSlopededx->SetMarkerStyle(22);
  grSlopededx->SetLineColor(2);
  grSlopededx->SetMarkerColor(2);
  grSlopededx->Draw("Psame");
  myc[0]->cd(15);
  grOffsetx0->SetTitle(";# layers removed;offset (mips)");
  grOffsetx0->SetMarkerStyle(20);
  grOffsetx0->Draw("AP");
  grOffsetdedx->SetMarkerStyle(22);
  grOffsetdedx->SetLineColor(2);
  grOffsetdedx->SetMarkerColor(2);
  grOffsetdedx->Draw("Psame");
  myc[0]->Update();
  myc[0]->Print(("Descop/Calibration_"+savestr+".pdf").c_str());

  myc[1]->cd(14);
  grResox0->SetTitle(";# layers removed;#sigma/E");
  grResox0->SetMarkerStyle(20);
  grResox0->SetMinimum(0.017);
  grResox0->SetMaximum(0.03);
  grResox0->Draw("AP");
  grResodedx->SetMarkerStyle(22);
  grResodedx->SetLineColor(2);
  grResodedx->SetMarkerColor(2);
  grResodedx->Draw("Psame");
  myc[1]->cd(15);
  grResoCalibx0->SetTitle(";# layers removed;#sigma/E");
  grResoCalibx0->SetMarkerStyle(20);
  grResoCalibx0->SetMinimum(0.017);
  grResoCalibx0->SetMaximum(0.03);
  grResoCalibx0->Draw("AP");
  grResoCalibdedx->SetMarkerStyle(22);
  grResoCalibdedx->SetLineColor(2);
  grResoCalibdedx->SetMarkerColor(2);
  grResoCalibdedx->Draw("Psame");
  myc[1]->Update();
  myc[1]->Print(("Descop/ResolutionRaw_"+savestr+".pdf").c_str());

  //if (doCMSSW){
  myc[2]->cd(14);
  grResox0->Draw("AP");
  grResodedx->Draw("Psame");
  myc[2]->cd(15);
  grResoCalibx0->Draw("AP");
  grResoCalibdedx->Draw("Psame");
  myc[2]->Update();
  myc[2]->Print(("Descop/ResolutionCalib_"+savestr+".pdf").c_str());
    //}

  //myc[3]->cd(1);
  //grReso->Draw("AP");
  myc[3]->cd();
  grResoCalibx0->Draw("AP");
  grResoCalibdedx->Draw("Psame");
  myc[3]->Update();
  myc[3]->Print(("Descop/ResolutionvsnRemove_"+savestr+".pdf").c_str());










  return 0;
}
