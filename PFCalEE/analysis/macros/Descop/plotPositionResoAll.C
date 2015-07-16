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

std::string plotDir(const bool doPhi45,
		    const unsigned version,
		    const unsigned etabin,
		    const unsigned pt,
		    const unsigned pu){
  std::ostringstream dir;
  if (!doPhi45) {
    if (version==12){
      if ((pu<200 && etabin!=19)) 
	dir << "../PLOTS/gitV00-02-12/version12/gamma/200um/eta" << etabin << "_et" << pt << "_pu" << pu;
      else 
	dir << "/afs/cern.ch/work/a/amagnan/PFCalEEAna/PLOTS/gitV00-02-12/version12/gamma/200um/eta" << etabin << "_et" << pt << "_pu" << pu;
    }
    else {
      dir << "/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalDescopPos/gitV05-02-04/version" << version << "/gamma/fixedVtx_/eta" << etabin << "_et" << pt << "_pu" << pu << "_IC3";
    }
  } else {
    dir << "/afs/cern.ch/work/a/amagnan/PFCalEEAna/PLOTS/gitV00-02-12/version12/gamma/200um/phi_0.250pi/eta" << etabin << "_et" << pt << "_pu" << pu;
  }
  return dir.str();
};

double E(const unsigned pT, const unsigned eta){
  return pT*cosh(eta/10.);
};

double pT(const unsigned E, const unsigned eta){
  return E/cosh(eta/10.);
};

struct Result{
  unsigned nEvts;
  double fitQuality;
  double fitQuality_rms;
  double impactX;
  double impactY;
  double tanAngleX;
  double tanAngleY;
  double impactX_rms;
  double impactY_rms;
  double tanAngleX_rms;
  double tanAngleY_rms;
  double residual_xFF;
  double residual_yFF;
  double residual_x14;
  double residual_y14;
  double residual_tanx;
  double residual_tany;
  double residual_xFFerr;
  double residual_yFFerr;
  double residual_x14err;
  double residual_y14err;
  double residual_tanxerr;
  double residual_tanyerr;
   double residual_xFF_rms;
  double residual_yFF_rms;
  double residual_x14_rms;
  double residual_y14_rms;
  double residual_tanx_rms;
  double residual_tany_rms;
  double residual_xFF_rmserr;
  double residual_yFF_rmserr;
  double residual_x14_rmserr;
  double residual_y14_rmserr;
  double residual_tanx_rmserr;
  double residual_tany_rmserr;
};

void fillPuPlots(std::string plotDir,unsigned etabin, unsigned pt, unsigned pu, TH1F *& meanPu, TH1F *& eventPu, TH1F *& diffPu){

  TFile *file = TFile::Open((plotDir+".root").c_str());
  
  if (!file){
    std::cout << " -- Error, input file " << plotDir << ".root cannot be opened. Skipping..." << std::endl;
    return; 
  }
  else std::cout << " -- file " << file->GetName() << " successfully opened." << std::endl;

  std::string suffix = "";
  
  file->cd("PositionFit");

  SetTdrStyle();
  gStyle->SetOptStat(0);

  TH2F *p_hitMeanPuContrib = (TH2F*)gDirectory->Get("p_hitMeanPuContrib");
  TH2F *p_hitEventPuContrib = (TH2F*)gDirectory->Get("p_hitEventPuContrib");
  TH1F *p_diffPuContrib = (TH1F*)gDirectory->Get("p_diffPuContrib");

  //SR2 calib
  double slope = 79.;

  if (p_hitMeanPuContrib && p_hitEventPuContrib && p_diffPuContrib){

    std::ostringstream lname;
    lname << "meanPu_eta" << etabin << "_pu" << pu;
    if (!meanPu) meanPu = (TH1F*)p_hitMeanPuContrib->ProjectionY(lname.str().c_str());
    else meanPu->Add(p_hitMeanPuContrib->ProjectionY());
    std::cout << " -- meanPu entries: " << meanPu->GetEntries() << std::endl;
    lname.str("");
    lname << "eventPu_eta" << etabin << "_pu" << pu;
    if (!eventPu) eventPu = (TH1F*)p_hitEventPuContrib->ProjectionY(lname.str().c_str());
    else eventPu->Add(p_hitEventPuContrib->ProjectionY());
    std::cout << " -- eventPu entries: " << eventPu->GetEntries() << std::endl;

    for (int bin(1);bin < p_diffPuContrib->GetNbinsX()+1;++bin){
      diffPu->Fill(p_diffPuContrib->GetBinCenter(bin)/slope,p_diffPuContrib->GetBinContent(bin));
    }
    std::cout << " -- diffPu entries: " << diffPu->Integral() << std::endl;
  }
  
}

Result plotOnePoint(std::string plotDir, TCanvas * & mycFit, TString sumDir, unsigned etabin, unsigned pt, unsigned pu){

  Result res;
  res.nEvts = 0;
  res.fitQuality = 0;
  res.fitQuality_rms = 0;
  res.impactX = 0;
  res.impactY = 0;
  res.tanAngleX = 0;
  res.tanAngleY = 0;
  res.impactX_rms = 0;
  res.impactY_rms = 0;
  res.tanAngleX_rms = 0;
  res.tanAngleY_rms = 0;
  res.residual_xFF = 0;
  res.residual_yFF = 0;
  res.residual_x14 = 0;
  res.residual_y14 = 0;
  res.residual_tanx = 0;
  res.residual_tany = 0;
  res.residual_xFFerr = 0;
  res.residual_yFFerr = 0;
  res.residual_x14err = 0;
  res.residual_y14err = 0;
  res.residual_tanxerr = 0;
  res.residual_tanyerr = 0;
  res.residual_xFF_rms = 0;
  res.residual_yFF_rms = 0;
  res.residual_x14_rms = 0;
  res.residual_y14_rms = 0;
  res.residual_tanx_rms = 0;
  res.residual_tany_rms = 0;
  res.residual_xFF_rmserr = 0;
  res.residual_yFF_rmserr = 0;
  res.residual_x14_rmserr = 0;
  res.residual_y14_rmserr = 0;
  res.residual_tanx_rmserr = 0;
  res.residual_tany_rmserr = 0;

  TFile *file = TFile::Open((plotDir+".root").c_str());
  
  if (!file){
    std::cout << " -- Error, input file " << plotDir << ".root cannot be opened. Skipping..." << std::endl;
    return res;
  }
  else std::cout << " -- file " << file->GetName() << " successfully opened." << std::endl;

  std::string suffix = "";
  
  file->cd("PositionFit");

  SetTdrStyle();
  gStyle->SetOptStat("eMRuo");

  //TCanvas *myc = new TCanvas("myc","myc",1);

  TH2F *p_hitMeanPuContrib = (TH2F*)gDirectory->Get("p_hitMeanPuContrib");
  TH2F *p_hitEventPuContrib = (TH2F*)gDirectory->Get("p_hitEventPuContrib");

  if (p_hitMeanPuContrib && p_hitEventPuContrib){

    TCanvas *mycE = new TCanvas("mycE","mycE",1500,750);//1500,1000);

    mycE->Divide(2,1);
    mycE->cd(1);
    gPad->SetRightMargin(0.15);
    gPad->SetLogz(1);
    p_hitMeanPuContrib->SetStats(0);
    p_hitMeanPuContrib->GetYaxis()->SetRangeUser(0,7);
    p_hitMeanPuContrib->GetZaxis()->SetTitleOffset(0.5);
    p_hitMeanPuContrib->Draw("colz");
    mycE->cd(2);
    gPad->SetRightMargin(0.15);
    gPad->SetLogz(1);
    p_hitEventPuContrib->SetStats(0);
    p_hitEventPuContrib->GetYaxis()->SetRangeUser(0,7);
    p_hitEventPuContrib->GetZaxis()->SetTitleOffset(0.5);
    p_hitEventPuContrib->Draw("colz");
    
    mycE->Update();
    mycE->Print((plotDir+"/AveragePuE"+suffix+".pdf").c_str());
  }
  //return 1;

  TCanvas *mycL = new TCanvas("mycL","mycL",1500,1000);
  TCanvas *mycD = new TCanvas("mycD","mycD",1500,1000);
  TCanvas *mycR = new TCanvas("mycR","mycR",1500,1000);
  TCanvas *mycW = new TCanvas("mycW","mycW",1);
  TLatex lat;
  char buf[500];

  TH2D *p_errorMatrix_x = (TH2D*)gDirectory->Get("p_errorMatrix_x");
  TH2D *p_corrMatrix_x = (TH2D*)gDirectory->Get("p_corrMatrix_x");
  TH2D *p_errorMatrix_y = (TH2D*)gDirectory->Get("p_errorMatrix_y");
  TH2D *p_corrMatrix_y = (TH2D*)gDirectory->Get("p_corrMatrix_y");
  bool skipDetailed = false;
  if (!p_errorMatrix_x) {
    std::cout << " -- Warning, input file " << plotDir << " does not contain detailed fit histos..." << std::endl;
    //return res;
    skipDetailed = true;
    return res;
  }

  TH1F *p_chi2overNDF = (TH1F*)gDirectory->Get("p_chi2overNDF");

  if (!p_chi2overNDF){
    std::cout << " -- Error, input file " << plotDir << ".root does not contain fit histos. Skipping..." << std::endl;
    return res;
  }


  TH1F *p_impactX = (TH1F*)gDirectory->Get("p_impactXFF");
  TH1F *p_impactXFF_residual = (TH1F*)gDirectory->Get("p_impactXFF_residual");
  TH1F *p_impactX14_residual = (TH1F*)gDirectory->Get("p_impactX14_residual");
  TH1F *p_impactX_truth = (TH1F*)gDirectory->Get("p_impactXFF_truth");
  TH1F *p_tanAngleX = (TH1F*)gDirectory->Get("p_tanAngleX");
  TH1F *p_tanAngleX_residual = (TH1F*)gDirectory->Get("p_tanAngleX_residual");
  TH1F *p_angleX_residual = (TH1F*)gDirectory->Get("p_angleX_residual");
  TH1F *p_tanAngleX_truth = (TH1F*)gDirectory->Get("p_tanAngleX_truth");
  TH1F *p_impactY = (TH1F*)gDirectory->Get("p_impactYFF");
  TH1F *p_impactYFF_residual = (TH1F*)gDirectory->Get("p_impactYFF_residual");
  TH1F *p_impactY14_residual = (TH1F*)gDirectory->Get("p_impactY14_residual");
  TH1F *p_impactY_truth = (TH1F*)gDirectory->Get("p_impactYFF_truth");
  TH1F *p_tanAngleY = (TH1F*)gDirectory->Get("p_tanAngleY");
  TH1F *p_tanAngleY_residual = (TH1F*)gDirectory->Get("p_tanAngleY_residual");
  TH1F *p_angleY_residual = (TH1F*)gDirectory->Get("p_angleY_residual");
  TH1F *p_tanAngleY_truth = (TH1F*)gDirectory->Get("p_tanAngleY_truth");
  
  TH1F *p_nLayersFit = (TH1F*)gDirectory->Get("p_nLayersFit");
  //TH2F *p_etavsphi_max = (TH2F*)gDirectory->Get("p_etavsphi_max");
  TH2F *p_recoXvsLayer = (TH2F*)gDirectory->Get("p_recoXvsLayer");
  TH2F *p_recoYvsLayer = (TH2F*)gDirectory->Get("p_recoYvsLayer");
  TH2F *p_recoZvsLayer = (TH2F*)gDirectory->Get("p_recoZvsLayer");
  TH2F *p_truthXvsLayer = (TH2F*)gDirectory->Get("p_truthXvsLayer");
  TH2F *p_truthYvsLayer = (TH2F*)gDirectory->Get("p_truthYvsLayer");
  TH2F *p_fitXvsLayer = (TH2F*)gDirectory->Get("p_fitXvsLayer");
  TH2F *p_fitYvsLayer = (TH2F*)gDirectory->Get("p_fitYvsLayer");
  TH1F *p_positionReso = (TH1F*)gDirectory->Get("p_positionReso");
  TH1F *p_angularReso = (TH1F*)gDirectory->Get("p_angularReso");

  mycL->Divide(4,2);
  mycL->cd(1);
  gPad->SetLogz(1);
  if (!skipDetailed){
    p_errorMatrix_x->SetStats(0);
    p_errorMatrix_x->SetMinimum(0.01);
    p_errorMatrix_x->Draw("colz");
  }
  mycL->cd(5);
  gPad->SetLogz(1);
  if (!skipDetailed){
    p_corrMatrix_x->SetStats(0);
    p_corrMatrix_x->SetMinimum(0.01);
    p_corrMatrix_x->Draw("colz");
  }
  mycL->cd(2);
  gPad->SetLogz(1);
  if (!skipDetailed){
    p_errorMatrix_y->SetStats(0);
    p_errorMatrix_y->SetMinimum(0.01);
    p_errorMatrix_y->Draw("colz");
  }
  mycL->cd(6);
  gPad->SetLogz(1);
  if (!skipDetailed){
    p_corrMatrix_y->SetStats(0);
    p_corrMatrix_y->SetMinimum(0.01);
    p_corrMatrix_y->Draw("colz");
  }

  mycL->cd(3);
  //gPad->SetLogy(1);
  p_impactX->SetLineColor(1);
  p_impactX->SetMarkerColor(1);
  p_impactX->SetMarkerStyle(21);
  p_impactX->StatOverflows(0);
  double minX = p_impactX_truth->GetMean()-20;
  double maxX = p_impactX_truth->GetMean()+20;
  p_impactX->GetXaxis()->SetRangeUser(minX,maxX);
  //p_impactX->SetMaximum(p_impactX_truth->GetMaximum()*1.1);
  p_impactX->Draw("PE");
  p_impactX_truth->SetLineColor(2);
  p_impactX_truth->Draw("same");

  res.nEvts = p_impactX->GetEntries();
  res.impactX = p_impactX->GetMean();
  res.impactX_rms = p_impactX->GetRMS();

  mycL->cd(4);
  //gPad->SetLogy(1);
  p_impactY->SetLineColor(1);
  p_impactY->SetMarkerColor(1);
  p_impactY->SetMarkerStyle(21);
  p_impactY->StatOverflows(0);
  double minY = p_impactY_truth->GetMean()-20;
  double maxY = p_impactY_truth->GetMean()+60;
  p_impactY->GetXaxis()->SetRangeUser(minY,maxY);
  //p_impactY->SetMaximum(p_impactY_truth->GetMaximum()*1.1);
  p_impactY->Draw("PE");
  p_impactY_truth->SetLineColor(2);
  p_impactY_truth->Draw("same");

  res.impactY = p_impactY->GetMean();
  res.impactY_rms = p_impactY->GetRMS();

  mycL->cd(7);
  gPad->SetLogy(1);
  p_tanAngleX->SetLineColor(1);
  p_tanAngleX->SetMarkerColor(1);
  p_tanAngleX->SetMarkerStyle(21);
  p_tanAngleX->GetXaxis()->SetRangeUser(p_tanAngleX_truth->GetMean()-0.2,p_tanAngleX_truth->GetMean()+0.2);
  p_tanAngleX->SetMaximum(p_tanAngleX_truth->GetMaximum()*1.1);
  p_tanAngleX->Draw("PE");
  p_tanAngleX_truth->SetLineColor(2);
  p_tanAngleX_truth->Draw("same");

  res.tanAngleX = p_tanAngleX->GetMean();
  res.tanAngleX_rms = p_tanAngleX->GetRMS();


  mycL->cd(8);
  gPad->SetLogy(1);
  p_tanAngleY->SetLineColor(1);
  p_tanAngleY->SetMarkerColor(1);
  p_tanAngleY->SetMarkerStyle(21);
  p_tanAngleY->GetXaxis()->SetRangeUser(p_tanAngleY_truth->GetMean()-0.2,p_tanAngleY_truth->GetMean()+0.2);
  p_tanAngleY->SetMaximum(p_tanAngleY_truth->GetMaximum()*1.1);
  p_tanAngleY->Draw("PE");
  p_tanAngleY_truth->SetLineColor(2);
  p_tanAngleY_truth->Draw("same");

  res.tanAngleY = p_tanAngleY->GetMean();
  res.tanAngleY_rms = p_tanAngleY->GetRMS();

  mycL->Update();
  mycL->Print((plotDir+"/PositionFitSummary"+suffix+".pdf").c_str());

  mycD->Divide(3,2);
  mycD->cd(1);
  gPad->SetLogy(1);
  p_nLayersFit->Draw();
  //p_etavsphi_max->Draw("colz");
  //p_etavsphi_max->SetStats(0);
  mycD->cd(4);
  gPad->SetLogy(1);
  gStyle->SetStatW(0.3);
  //gStyle->SetStatH(0.3);
  p_chi2overNDF->GetXaxis()->SetRangeUser(0,20);
  p_chi2overNDF->Draw();

  res.fitQuality = p_chi2overNDF->GetMean();
  res.fitQuality_rms = p_chi2overNDF->GetRMS();
  //p_recoZvsLayer->Draw("colz");
  //p_recoZvsLayer->SetStats(0);
  
  mycD->cd(2);
  gPad->SetLogz(1);
  p_recoXvsLayer->RebinY(2);
  p_truthXvsLayer->RebinY(2);
  p_recoXvsLayer->GetYaxis()->SetRangeUser(p_recoXvsLayer->GetMean(2)-100,p_recoXvsLayer->GetMean(2)+100);
  p_recoXvsLayer->Draw("colz");
  p_recoXvsLayer->SetStats(0);
  p_truthXvsLayer->SetMarkerStyle(1);
  p_truthXvsLayer->Draw("same");
  p_truthXvsLayer->SetStats(0);
  mycD->cd(5);
  gPad->SetLogz(1);
  p_recoYvsLayer->RebinY(4);
  p_truthYvsLayer->RebinY(4);
  p_recoYvsLayer->GetYaxis()->SetRangeUser(p_recoYvsLayer->GetMean(2)-200,p_recoYvsLayer->GetMean(2)+200);
  p_recoYvsLayer->Draw("colz");
  p_recoYvsLayer->SetStats(0);
  p_truthYvsLayer->SetMarkerStyle(1);
  p_truthYvsLayer->Draw("same");
  p_truthYvsLayer->SetStats(0);
  mycD->cd(3);
  gPad->SetLogz(1);
  //p_fitXvsLayer->RebinY(2);
  p_fitXvsLayer->GetYaxis()->SetRangeUser(p_fitXvsLayer->GetMean(2)-20,p_fitXvsLayer->GetMean(2)+20);
  p_fitXvsLayer->Draw("colz");
  p_fitXvsLayer->SetStats(0);
  mycD->cd(6);
  gPad->SetLogz(1);
  //p_fitYvsLayer->RebinY(4);
  p_fitYvsLayer->GetYaxis()->SetRangeUser(p_fitYvsLayer->GetMean(2)-100,p_fitYvsLayer->GetMean(2)+100);
  p_fitYvsLayer->Draw("colz");
  p_fitYvsLayer->SetStats(0);

  mycD->Update();
  mycD->Print((plotDir+"/PositionFitDebug"+suffix+".pdf").c_str());

  mycR->Divide(4,2);
  mycR->cd(1);
  p_positionReso->GetXaxis()->SetRangeUser(p_positionReso->GetMean()-5*p_positionReso->GetRMS(),
					    p_positionReso->GetMean()+5*p_positionReso->GetRMS());
  p_positionReso->Draw();
  mycR->cd(5);
  p_angularReso->GetXaxis()->SetRangeUser(0,0.1);
  //p_angularReso->GetMean()-5*p_angularReso->GetRMS(),
  //p_angularReso->GetMean()+5*p_angularReso->GetRMS());
  p_angularReso->Draw();
  mycR->cd(2);
  p_impactXFF_residual->GetXaxis()->SetRangeUser(-2.5,5);
  p_impactXFF_residual->Draw();
  mycR->cd(3);
  p_impactX14_residual->GetXaxis()->SetRangeUser(-2.5,5);
  p_impactX14_residual->Draw();

  mycFit->cd();
  p_impactXFF_residual->Draw("PE");
  p_impactXFF_residual->Fit("gaus","0+");
  TF1 *fit = p_impactXFF_residual->GetFunction("gaus");
  if (!fit) return res;
  fit->SetLineColor(6);
  fit->Draw("same");
  res.residual_xFF = fit->GetParameter(1);
  res.residual_xFFerr = fit->GetParError(1);
  res.residual_xFF_rms = fit->GetParameter(2);
  res.residual_xFF_rmserr = fit->GetParError(2);
  sprintf(buf,"#eta=%3.1f, pT=%d GeV, pu=%d",etabin/10.,pt,pu);
  lat.DrawLatexNDC(0.1,0.96,buf);
  mycFit->Print(sumDir+"/Positionfits_x.pdf");

  if (etabin==19 && pt==50 && pu==0) mycFit->Print(sumDir+"/posFF_x_eta19_et50_pu0.C");
  if (etabin==19 && pt==50 && pu==140) mycFit->Print(sumDir+"/posFF_x_eta19_et50_pu140.C");

  mycFit->cd();
  p_impactX14_residual->Draw("PE");
  p_impactX14_residual->Fit("gaus","0+");
  fit = p_impactX14_residual->GetFunction("gaus");
  if (!fit) return res;
  fit->SetLineColor(6);
  fit->Draw("same");
  res.residual_x14 = fit->GetParameter(1);
  res.residual_x14err = fit->GetParError(1);
  res.residual_x14_rms = fit->GetParameter(2);
  res.residual_x14_rmserr = fit->GetParError(2);
  sprintf(buf,"#eta=%3.1f, pT=%d GeV, pu=%d",etabin/10.,pt,pu);
  lat.DrawLatexNDC(0.1,0.96,buf);
  mycFit->Print(sumDir+"/Positionfits_x.pdf");
  if (etabin==19 && pt==50 && pu==0) mycFit->Print(sumDir+"/pos14_x_eta19_et50_pu0.C");
  if (etabin==19 && pt==50 && pu==140) mycFit->Print(sumDir+"/pos14_x_eta19_et50_pu140.C");

  mycR->cd(6);
  p_impactYFF_residual->GetXaxis()->SetRangeUser(-2.5,5);
  p_impactYFF_residual->Draw();
  mycR->cd(7);
  p_impactY14_residual->GetXaxis()->SetRangeUser(-2.5,5);
  p_impactY14_residual->Draw();

  mycFit->cd();
  p_impactYFF_residual->Draw("PE");
  p_impactYFF_residual->Fit("gaus","0+");
  fit = p_impactYFF_residual->GetFunction("gaus");
  if (!fit) return res;
  fit->SetLineColor(6);
  fit->Draw("same");
  res.residual_yFF = fit->GetParameter(1);
  res.residual_yFFerr = fit->GetParError(1);
  res.residual_yFF_rms = fit->GetParameter(2);
  res.residual_yFF_rmserr = fit->GetParError(2);
  sprintf(buf,"#eta=%3.1f, pT=%d GeV, pu=%d",etabin/10.,pt,pu);
  lat.DrawLatexNDC(0.1,0.96,buf);
  mycFit->Print(sumDir+"/Positionfits_y.pdf");
  if (etabin==19 && pt==50 && pu==0) mycFit->Print(sumDir+"/posFF_y_eta19_et50_pu0.C");
  if (etabin==19 && pt==50 && pu==140) mycFit->Print(sumDir+"/posFF_y_eta19_et50_pu140.C");

  mycFit->cd();
  p_impactY14_residual->Draw("PE");
  p_impactY14_residual->Fit("gaus","0+");
  fit = p_impactY14_residual->GetFunction("gaus");
  if (!fit) return res;
  fit->SetLineColor(6);
  fit->Draw("same");
  res.residual_y14 = fit->GetParameter(1);
  res.residual_y14err = fit->GetParError(1);
  res.residual_y14_rms = fit->GetParameter(2);
  res.residual_y14_rmserr = fit->GetParError(2);
  sprintf(buf,"#eta=%3.1f, pT=%d GeV, pu=%d",etabin/10.,pt,pu);
  lat.DrawLatexNDC(0.1,0.96,buf);
  mycFit->Print(sumDir+"/Positionfits_y.pdf");
  if (etabin==19 && pt==50 && pu==0) mycFit->Print(sumDir+"/pos14_y_eta19_et50_pu0.C");
  if (etabin==19 && pt==50 && pu==140) mycFit->Print(sumDir+"/pos14_y_eta19_et50_pu140.C");

  mycR->cd(4);
  p_tanAngleX_residual->GetXaxis()->SetRangeUser(-0.03,0.06);
  p_tanAngleX_residual->Draw();

  mycFit->cd();
  p_angleX_residual->Draw("PE");
  p_angleX_residual->Fit("gaus","0+");
  fit = p_angleX_residual->GetFunction("gaus");
  if (!fit) return res;
  fit->SetLineColor(6);
  fit->Draw("same");
  res.residual_tanx = fit->GetParameter(1);
  res.residual_tanxerr = fit->GetParError(1);
  res.residual_tanx_rms = fit->GetParameter(2);//p_tanAngleX->GetRMS();
  res.residual_tanx_rmserr = fit->GetParError(2);//p_tanAngleX->GetRMS();
  sprintf(buf,"#eta=%3.1f, pT=%d GeV, pu=%d",etabin/10.,pt,pu);
  lat.DrawLatexNDC(0.1,0.96,buf);
  mycFit->Print(sumDir+"/Positionfits_x.pdf");
  if (etabin==19 && pt==50 && pu==0) mycFit->Print(sumDir+"/angle_x_eta19_et50_pu0.C");
  if (etabin==19 && pt==50 && pu==140) mycFit->Print(sumDir+"/angle_x_eta19_et50_pu140.C");



  mycR->cd(8);
  p_angleY_residual->GetXaxis()->SetRangeUser(-0.03,0.06);
  p_angleY_residual->Draw();
  mycFit->cd();
  p_angleY_residual->Draw("PE");
  p_angleY_residual->Fit("gaus","0+");
  if (!fit) return res;
  fit = p_angleY_residual->GetFunction("gaus");
  fit->SetLineColor(6);
  fit->Draw("same");
  res.residual_tany = fit->GetParameter(1);
  res.residual_tanyerr = fit->GetParError(1);
  res.residual_tany_rms = fit->GetParameter(2);//p_tanAngleX->GetRMS();
  res.residual_tany_rmserr = fit->GetParError(2);//p_tanAngleX->GetRMS();
  sprintf(buf,"#eta=%3.1f, pT=%d GeV, pu=%d",etabin/10.,pt,pu);
  lat.DrawLatexNDC(0.1,0.96,buf);
  mycFit->Print(sumDir+"/Positionfits_y.pdf");
  if (etabin==19 && pt==50 && pu==0) mycFit->Print(sumDir+"/angle_y_eta19_et50_pu0.C");
  if (etabin==19 && pt==50 && pu==140) mycFit->Print(sumDir+"/angle_y_eta19_et50_pu140.C");

  //res.residual_x = p_impactX_residual->GetMean();
  //res.residual_x_rms = p_impactX_residual->GetRMS();
  //res.residual_y = p_impactY_residual->GetMean();
  //res.residual_y_rms = p_impactY_residual->GetRMS();

  //res.residual_tanx = p_tanAngleX_residual->GetMean();
  //res.residual_tanx_rms = p_tanAngleX_residual->GetRMS();
  //res.residual_tany = p_tanAngleY_residual->GetMean();
  //res.residual_tany_rms = p_tanAngleY_residual->GetRMS();


  mycR->Update();
  mycR->Print((plotDir+"/PositionFitReso"+suffix+".pdf").c_str());

  return res;

}//main

TPad* plot_ratio(TCanvas *canv, unsigned nb){
  canv->SetFillColor      (0);
  canv->SetBorderMode     (0);
  canv->SetBorderSize     (10);
  // Set margins to reasonable defaults
  canv->SetLeftMargin     (0.18);
  canv->SetLeftMargin     (0.17);
  canv->SetRightMargin    (0.05);
  canv->SetTopMargin      (0.12);
  canv->SetBottomMargin   (0.18);
  // Setup a frame which makes sense
  canv->SetFrameFillStyle (0);
  canv->SetFrameLineStyle (0);
  canv->SetFrameBorderMode(0);
  canv->SetFrameBorderSize(10);
  canv->SetFrameFillStyle (0);
  canv->SetFrameLineStyle (0);
  canv->SetFrameBorderMode(0);
  canv->SetFrameBorderSize(10);      

  canv->cd();
  TPad *pad = 0;
  const unsigned neta=7;
  std::ostringstream padname;
  padname << "eta"<<nb;
  pad = new TPad(padname.str().c_str(),"pad",0, (neta-nb-1)*1./neta ,1 ,(neta-nb)*1./neta);
  if (nb==0) {
    pad->SetBottomMargin(0.05);
    pad->SetTopMargin(0.09);
  }
  else if (nb==neta-1){
    pad->SetTopMargin(0.05);
    pad->SetBottomMargin(0.24);
  }
  else {
    pad->SetTopMargin(0.05);
    pad->SetBottomMargin(0.05);
  }
  pad->Draw();
  pad->cd();
  return pad;
  
};

int plotPositionResoAll(){//main

  SetTdrStyle();

  const bool doPhi45 = false;
  const unsigned nV = 3;
  const unsigned version[nV] = {30,34,35};

  const unsigned npu = 1;//3;
  const unsigned neta = 6;//4;//(doPhi45 || version[iv]==13) ? 4 : 7;

  const bool doVsE = true;

  unsigned pu[npu] = {0};//,200};
  unsigned etaPatch[neta] = {17,19,21,23,25,27};
  //unsigned etaPatch[neta] = {17,21,23,25,27,29};
  unsigned eta[neta] = {17,19,21,23,25,27};
  //for (unsigned ieta(0);ieta<neta;++ieta){
  //eta[ieta] = 17+4*ieta;
    //if (doPhi45 || version[iv]==13) eta[ieta] = 17+4*ieta;
    //else eta[ieta] = 17+2*ieta;
  //}
  const unsigned npt = 1;
  unsigned pt[npt] = {70};//{3,5,7,10,20,30,40,50,60,70,80,90,100,125,150,175,200};


  TCanvas *mycFit = new TCanvas("mycFit","mycFit",1);
  

  TCanvas *mycNice = new TCanvas("mycNice","mycNice",1500,750);

  TCanvas *mycPu = new TCanvas("mycPu","mycPu",1500,1000);
  TCanvas *mycFitQual = new TCanvas("mycFitQual","mycFitQual",1500,1000);
  TCanvas *mycEvt = new TCanvas("mycEvt","mycEvt",1500,1000);
  TCanvas *mycRP14 = new TCanvas("mycRP14","mycRP14",1500,1000);
  TCanvas *mycRPFF = new TCanvas("mycRPFF","mycRPFF",1500,1000);
  TCanvas *mycRA = new TCanvas("mycRA","mycRA",1500,1000);
  TCanvas *mycRP14sig = new TCanvas("mycRP14sig","mycRP14sig",1500,1000);
  TCanvas *mycRPFFsig = new TCanvas("mycRPFFsig","mycRPFFsig",1500,1000);
  TCanvas *mycRAsig = new TCanvas("mycRAsig","mycRAsig",1500,1000);
  mycPu->Divide(2,1);
  //mycEvt->Divide(2,1);
  //mycFitQual->Divide(2,1);
  //mycRP14->Divide(2,1);
  //mycRPFF->Divide(2,1);
  //mycRA->Divide(2,1);
  // mycRP14sig->Divide(2,1);
  //mycRPFFsig->Divide(2,1);
  //mycRAsig->Divide(2,1);
  /*const unsigned nCan = 9;
  TPad *mypad[nCan][neta];
  TPad *left[nCan];
  TPad *right[nCan];
  left[0] = (TPad*)mycEvt->cd(1);
  right[0] = (TPad*)mycEvt->cd(2);
  left[1] = (TPad*)mycFitQual->cd(1);
  right[1] = (TPad*)mycFitQual->cd(2);
  left[2] = (TPad*)mycRP14->cd(1);
  right[2] = (TPad*)mycRP14->cd(2);
  left[3] = (TPad*)mycRA->cd(1);
  right[3] = (TPad*)mycRA->cd(2);
  left[5] = (TPad*)mycRPFF->cd(1);
  right[5] = (TPad*)mycRPFF->cd(2);
  left[4] = (TPad*)mycPu->cd(1);
  right[4] = (TPad*)mycPu->cd(2);
  left[6] = (TPad*)mycRPFFsig->cd(1);
  right[6] = (TPad*)mycRPFFsig->cd(2);
  left[7] = (TPad*)mycRP14sig->cd(1);
  right[7] = (TPad*)mycRP14sig->cd(2);
  left[8] = (TPad*)mycRAsig->cd(1);
  right[8] = (TPad*)mycRAsig->cd(2);
  for (unsigned iC(0);iC<nCan;++iC){
    left[iC]->Divide(1,3);
    right[iC]->Divide(1,3);
    //left[iC]->Divide(1,2);
    //right[iC]->Divide(1,2);
    for (unsigned ieta=0; ieta<neta;++ieta){//loop on pt values
      if (ieta< 3) mypad[iC][ieta] = (TPad*)left[iC]->GetPad(ieta+1);
      else mypad[iC][ieta] = (TPad*)right[iC]->GetPad((ieta-3)+1);
      //if (ieta< 2) mypad[iC][ieta] = (TPad*)left[iC]->GetPad(ieta+1);
      //else mypad[iC][ieta] = (TPad*)right[iC]->GetPad((ieta-2)+1);
    }
  }
  */

  TGraphErrors *grEvts[nV][npu];
  TGraphErrors *grFit[nV][npu];
  TGraphErrors *grResidualXFF[nV][npu];
  TGraphErrors *grResidualYFF[nV][npu];
  TGraphErrors *grResidualX14[nV][npu];
  TGraphErrors *grResidualY14[nV][npu];
  TGraphErrors *grResidualTanX[nV][npu];
  TGraphErrors *grResidualTanY[nV][npu];
  TGraphErrors *grResolXFF[nV][npu];
  TGraphErrors *grResolYFF[nV][npu];
  TGraphErrors *grResolX14[nV][npu];
  TGraphErrors *grResolY14[nV][npu];
  TGraphErrors *grResolTanX[nV][npu];
  TGraphErrors *grResolTanY[nV][npu];

  for (unsigned iv(0); iv<nV;++iv){//loop on versions

    TString sumDir = version[iv]==13? "PLOTS_v13" : doPhi45 ? "PLOTS_phi45" : "PLOTS_v" ;
    sumDir+= version[iv];
    sumDir += "_pT";
    sumDir += pt[0];
    mycFit->Print(sumDir+"/Positionfits_x.pdf[");
    mycFit->Print(sumDir+"/Positionfits_y.pdf[");

    double ptval[npu][npt];
    double etaval[npu][neta];
    //plot pu
    
    for (unsigned ipu(1); ipu<npu;++ipu){//loop on pu values
      TH1F *meanPuProj[neta];
      TH1F *eventPuProj[neta];
      TH1F *diffPu[neta];
      for (unsigned ieta=0; ieta<neta;++ieta){//loop on eta values
	meanPuProj[ieta] = 0;
	eventPuProj[ieta] = 0;
	std::ostringstream lname;
	lname.str("");
	lname << "diffPu_eta" << ieta << "_pu" << pu[ipu];
	diffPu[ieta] = new TH1F(lname.str().c_str(),";E_{PU}^{RC}-E_{PU}^{Mean} (GeV);events",100,-5,25);
	
	for (unsigned ipt(0);ipt<npt;++ipt){
	  fillPuPlots(plotDir(doPhi45,version[iv],etaPatch[ieta],pt[ipt],pu[ipu]),eta[ieta],pt[ipt],pu[ipu],meanPuProj[ieta],eventPuProj[ieta],diffPu[ieta]);
	}
	if (!meanPuProj[ieta] || !eventPuProj[ieta]) {
	  std::cout << " Pu histos not found, continuing..." << std::endl;
	  continue;
	}
	mycPu->cd();
	gStyle->SetOptStat(0);
	meanPuProj[ieta]->SetMarkerStyle(21);
	meanPuProj[ieta]->GetXaxis()->SetRangeUser(0,20);
	meanPuProj[ieta]->GetXaxis()->SetTitle("E^{hit}_{PU} (MIPs)");
	meanPuProj[ieta]->Draw("PE");
	eventPuProj[ieta]->SetLineColor(2);
	eventPuProj[ieta]->SetMarkerColor(2);
	eventPuProj[ieta]->SetMarkerStyle(22);
	eventPuProj[ieta]->Draw("PEsame");
	
	TLegend *leg3 = new TLegend(0.65,0.74,0.94,0.94);
	leg3->SetFillColor(10);
	leg3->AddEntry(meanPuProj[ieta],"Avg from expo fit","P");
	leg3->AddEntry(eventPuProj[ieta],"Evt-by-evt from RC","P");
	leg3->Draw("same");
	TLatex lat;
	lat.SetTextSize(0.08);
	char buf[500];
	sprintf(buf,"<%d> pu, #eta = %3.1f",pu[ipu],eta[ieta]/10.);
	lat.DrawLatexNDC(0.2,0.85,buf);
	sprintf(buf,"RMS <pu> = %3.1f",meanPuProj[ieta]->GetRMS());
	lat.DrawLatexNDC(0.2,0.65,buf);
	sprintf(buf,"RMS pu^{RC} = %3.1f",eventPuProj[ieta]->GetRMS());
	lat.DrawLatexNDC(0.2,0.5,buf);
	
	lat.DrawLatexNDC(0.02,0.02,"HGCAL Geant4 Standalone");
      }
      std::ostringstream lPrint;
      lPrint << sumDir << "/SummaryAll_Pu" << pu[ipu];
      lPrint << ".pdf";
      mycPu->Update();
      mycPu->Print(lPrint.str().c_str());
      
      for (unsigned ieta=0; ieta<neta;++ieta){//loop on eta values
	if (!diffPu[ieta]){
	  std::cout << " Pu histos not found, continuing..." << std::endl;
	  continue;
	}

	mycPu->cd();
	gStyle->SetOptStat("eMRuo");
	diffPu[ieta]->Draw();
	TLatex lat;
	lat.SetTextSize(0.08);
	char buf[500];
	sprintf(buf,"<%d> pu, #eta = %3.1f",pu[ipu],eta[ieta]/10.);
	lat.DrawLatexNDC(0.2,0.85,buf);
	lat.DrawLatexNDC(0.02,0.02,"HGCAL Geant4 Standalone");
      }
      lPrint.str("");
      lPrint << sumDir << "/SummaryAll_DiffPu" << pu[ipu];
      lPrint << ".pdf";
      mycPu->Update();
      mycPu->Print(lPrint.str().c_str());
    }//loop on pu val
  
    //return 1;
    for (unsigned ipt(0);ipt<npt;++ipt){
      for (unsigned ipu(0); ipu<npu;++ipu){//loop on pu values
      
	std::vector<Result> resVec;
      
	unsigned newp = 0;
	for (unsigned ieta=0; ieta<neta;++ieta){//loop on eta values
	  Result lres = plotOnePoint(plotDir(doPhi45,version[iv],etaPatch[ieta]+((pt[0]==20 && iv==0 && ieta>0)?2:0),pt[ipt],pu[ipu]),mycFit,sumDir,eta[ieta],pt[ipt],pu[ipu]);
	  if (lres.nEvts > 0){
	    etaval[ipu][newp] = eta[ieta]/10.;//+0.1*(2.*ipu-1.);
	    if (doVsE) etaval[ipu][newp] = E(pt[ipt],eta[ieta]);//+2.*(2.*ipu-1.);
	    //std::cout << pu[ipu] << " " << eta[ie] << " " << etaval[ipu][newp] << std::endl;
	    newp++;
	    resVec.push_back(lres);
	  }
	}
      
	const unsigned nP = resVec.size();
	double pterr[nP];
	double nEvts[nP];
	double fitQuality[nP];
	double fitQuality_rms[nP];
	double impactX[nP];
	double impactY[nP];
	double tanAngleX[nP];
	double tanAngleY[nP];
	double impactX_rms[nP];
	double impactY_rms[nP];
	double tanAngleX_rms[nP];
	double tanAngleY_rms[nP];
	double residual_xFF[nP];
	double residual_yFF[nP];
	double residual_x14[nP];
	double residual_y14[nP];
	double residual_tanx[nP];
	double residual_tany[nP];
	double residual_xFFerr[nP];
	double residual_yFFerr[nP];
	double residual_x14err[nP];
	double residual_y14err[nP];
	double residual_tanxerr[nP];
	double residual_tanyerr[nP];
	double residual_xFF_rms[nP];
	double residual_yFF_rms[nP];
	double residual_x14_rms[nP];
	double residual_y14_rms[nP];
	double residual_tanx_rms[nP];
	double residual_tany_rms[nP];
	double residual_xFF_rmserr[nP];
	double residual_yFF_rmserr[nP];
	double residual_x14_rmserr[nP];
	double residual_y14_rmserr[nP];
	double residual_tanx_rmserr[nP];
	double residual_tany_rmserr[nP];

	double max_angle = 0;
	double max_pos = 0;
	for (unsigned iP(0); iP<nP;++iP){
	  pterr[iP] = 0;
	  nEvts[iP] = resVec[iP].nEvts;
	  fitQuality[iP] = resVec[iP].fitQuality;
	  fitQuality_rms[iP] = resVec[iP].fitQuality_rms;
	  impactX[iP] = resVec[iP].impactX;
	  impactY[iP] = resVec[iP].impactY;
	  tanAngleX[iP] = resVec[iP].tanAngleX;
	  tanAngleY[iP] = resVec[iP].tanAngleY;
	  impactX_rms[iP] = resVec[iP].impactX_rms;
	  impactY_rms[iP] = resVec[iP].impactY_rms;
	  tanAngleX_rms[iP] = resVec[iP].tanAngleX_rms;
	  tanAngleY_rms[iP] = resVec[iP].tanAngleY_rms;
	  residual_xFF[iP] = resVec[iP].residual_xFF;
	  residual_yFF[iP] = resVec[iP].residual_yFF;
	  residual_x14[iP] = resVec[iP].residual_x14;
	  residual_y14[iP] = resVec[iP].residual_y14;
	  residual_tanx[iP] = resVec[iP].residual_tanx*1000;
	  residual_tany[iP] = resVec[iP].residual_tany*1000;
	  residual_xFFerr[iP] = resVec[iP].residual_xFFerr;
	  residual_yFFerr[iP] = resVec[iP].residual_yFFerr;
	  residual_x14err[iP] = resVec[iP].residual_x14err;
	  residual_y14err[iP] = resVec[iP].residual_y14err;
	  residual_tanxerr[iP] = resVec[iP].residual_tanxerr*1000;
	  residual_tanyerr[iP] = resVec[iP].residual_tanyerr*1000;
	  residual_xFF_rms[iP] = resVec[iP].residual_xFF_rms;
	  residual_yFF_rms[iP] = resVec[iP].residual_yFF_rms;
	  residual_x14_rms[iP] = resVec[iP].residual_x14_rms;
	  residual_y14_rms[iP] = resVec[iP].residual_y14_rms;
	  residual_tanx_rms[iP] = resVec[iP].residual_tanx_rms*1000;
	  residual_tany_rms[iP] = resVec[iP].residual_tany_rms*1000;
	  residual_xFF_rmserr[iP] = resVec[iP].residual_xFF_rmserr;
	  residual_yFF_rmserr[iP] = resVec[iP].residual_yFF_rmserr;
	  residual_x14_rmserr[iP] = resVec[iP].residual_x14_rmserr;
	  residual_y14_rmserr[iP] = resVec[iP].residual_y14_rmserr;
	  residual_tanx_rmserr[iP] = resVec[iP].residual_tanx_rmserr*1000;
	  residual_tany_rmserr[iP] = resVec[iP].residual_tany_rmserr*1000;
	  if ( (fabs(residual_tanx[iP])+residual_tanx_rms[iP])> max_angle) max_angle = fabs(residual_tanx[iP])+residual_tanx_rms[iP];
	  if ( (fabs(residual_tany[iP])+residual_tany_rms[iP])> max_angle) max_angle = fabs(residual_tany[iP])+residual_tany_rms[iP];
	  if ( (fabs(residual_xFF[iP])+residual_xFF_rms[iP])> max_pos) max_pos = fabs(residual_xFF[iP])+residual_xFF_rms[iP];
	  if ( (fabs(residual_yFF[iP])+residual_yFF_rms[iP])> max_pos) max_pos = fabs(residual_yFF[iP])+residual_yFF_rms[iP];
	}
      
	grEvts[iv][ipu] = new TGraphErrors(nP,etaval[ipu],nEvts,pterr,pterr);
	grFit[iv][ipu] = new TGraphErrors(nP,etaval[ipu],fitQuality,pterr,fitQuality_rms);

	grResidualXFF[iv][ipu] = new TGraphErrors(nP,etaval[ipu],residual_xFF,pterr,residual_xFFerr);
	grResidualYFF[iv][ipu] = new TGraphErrors(nP,etaval[ipu],residual_yFF,pterr,residual_yFFerr);
	grResidualX14[iv][ipu] = new TGraphErrors(nP,etaval[ipu],residual_x14,pterr,residual_x14err);
	grResidualY14[iv][ipu] = new TGraphErrors(nP,etaval[ipu],residual_y14,pterr,residual_y14err);
	grResidualTanX[iv][ipu] = new TGraphErrors(nP,etaval[ipu],residual_tanx,pterr,residual_tanxerr);
	grResidualTanY[iv][ipu] = new TGraphErrors(nP,etaval[ipu],residual_tany,pterr,residual_tanyerr);

	grResolXFF[iv][ipu] = new TGraphErrors(nP,etaval[ipu],residual_xFF_rms,pterr,residual_xFF_rmserr);
	grResolYFF[iv][ipu] = new TGraphErrors(nP,etaval[ipu],residual_yFF_rms,pterr,residual_yFF_rmserr);
	grResolX14[iv][ipu] = new TGraphErrors(nP,etaval[ipu],residual_x14_rms,pterr,residual_x14_rmserr);
	grResolY14[iv][ipu] = new TGraphErrors(nP,etaval[ipu],residual_y14_rms,pterr,residual_y14_rmserr);
	grResolTanX[iv][ipu] = new TGraphErrors(nP,etaval[ipu],residual_tanx_rms,pterr,residual_tanx_rmserr);
	grResolTanY[iv][ipu] = new TGraphErrors(nP,etaval[ipu],residual_tany_rms,pterr,residual_tany_rmserr);

	mycEvt->cd();
	gPad->SetGridy(1);
	grEvts[iv][ipu]->SetTitle(";#eta;N_{events}");
	if (doVsE) grEvts[iv][ipu]->SetTitle(";E (GeV);N_{events}");
	grEvts[iv][ipu]->SetLineColor(ipu+1);
	grEvts[iv][ipu]->SetMarkerColor(ipu+1);
	grEvts[iv][ipu]->SetMarkerStyle(ipu+21);
	grEvts[iv][ipu]->SetMinimum(0);
	grEvts[iv][ipu]->SetMaximum(5200);
	grEvts[iv][ipu]->Draw(ipu==0?"APE":"PE");
      
	mycFitQual->cd();
	gPad->SetGridy(1);
	grFit[iv][ipu]->SetTitle(";#eta;#chi^{2}/NDF");
	if (doVsE) grFit[iv][ipu]->SetTitle(";E (GeV);#chi^{2}/NDF");
	grFit[iv][ipu]->SetLineColor(ipu+1);
	grFit[iv][ipu]->SetMarkerColor(ipu+1);
	grFit[iv][ipu]->SetMarkerStyle(ipu+21);
	grFit[iv][ipu]->SetMinimum(0);
	grFit[iv][ipu]->SetMaximum(10);
	grFit[iv][ipu]->Draw(ipu==0?"APE":"PE");
      
	mycRP14->cd();
	gPad->SetGridy(1);
	grResidualX14[iv][ipu]->SetTitle(";#eta;x14_{reco}-x14_{truth} (mm)");
	if (doVsE) grResidualX14[iv][ipu]->SetTitle(";E (GeV);x14_{reco}-x14_{truth} (mm)");
	grResidualX14[iv][ipu]->SetLineColor(ipu+1);
	grResidualX14[iv][ipu]->SetMarkerColor(ipu+1);
	grResidualX14[iv][ipu]->SetMarkerStyle(ipu+21);
	//grResidualX14[iv][ipu]->SetMaximum(4);//max_pos);
	//grResidualX14[iv][ipu]->SetMinimum(-4);//max_pos);
	grResidualX14[iv][ipu]->Draw(ipu==0?"APE":"PE");
      
	grResidualY14[iv][ipu]->SetLineColor(ipu+3);
	grResidualY14[iv][ipu]->SetMarkerColor(ipu+3);
	grResidualY14[iv][ipu]->SetMarkerStyle(ipu+23);
	grResidualY14[iv][ipu]->Draw("PE");
      
	mycRA->cd();
	gPad->SetGridy(1);
	grResidualTanX[iv][ipu]->SetTitle(";#eta;tanA_{reco}-tanA_{truth} (mrad)");
	if (doVsE) grResidualTanX[iv][ipu]->SetTitle(";E (GeV);tanA_{reco}-tanA_{truth} (mrad)");
	grResidualTanX[iv][ipu]->SetLineColor(ipu+1);
	grResidualTanX[iv][ipu]->SetMarkerColor(ipu+1);
	grResidualTanX[iv][ipu]->SetMarkerStyle(ipu+21);
	//grResidualTanX[iv][ipu]->SetMaximum(40);//max_angle);
	//grResidualTanX[iv][ipu]->SetMinimum(-40);//max_angle);
	grResidualTanX[iv][ipu]->Draw(ipu==0?"APE":"PE");
      
	grResidualTanY[iv][ipu]->SetLineColor(ipu+3);
	grResidualTanY[iv][ipu]->SetMarkerColor(ipu+3);
	grResidualTanY[iv][ipu]->SetMarkerStyle(ipu+23);
	grResidualTanY[iv][ipu]->Draw("PE");
      
	mycRPFF->cd();
	gPad->SetGridy(1);
	grResidualXFF[iv][ipu]->SetTitle(";#eta;xFF_{reco}-xFF_{truth} (mm)");
	if (doVsE) grResidualXFF[iv][ipu]->SetTitle(";E (GeV);xFF_{reco}-xFF_{truth} (mm)");
	grResidualXFF[iv][ipu]->SetLineColor(ipu+1);
	grResidualXFF[iv][ipu]->SetMarkerColor(ipu+1);
	grResidualXFF[iv][ipu]->SetMarkerStyle(ipu+21);
	//grResidualXFF[iv][ipu]->SetMaximum(4);//max_pos);
	//grResidualXFF[iv][ipu]->SetMinimum(-4);//max_pos);
	grResidualXFF[iv][ipu]->Draw(ipu==0?"APE":"PE");
      
	grResidualYFF[iv][ipu]->SetLineColor(ipu+3);
	grResidualYFF[iv][ipu]->SetMarkerColor(ipu+3);
	grResidualYFF[iv][ipu]->SetMarkerStyle(ipu+23);
	grResidualYFF[iv][ipu]->Draw("PE");
      
	mycRPFFsig->cd();
	gPad->SetGridy(1);
	//gPad->SetLogy(1);
	grResolXFF[iv][ipu]->SetTitle(";#eta;#sigma(xFF_{reco}-xFF_{truth}) (mm)");
	if (doVsE) grResolXFF[iv][ipu]->SetTitle(";E (GeV);#sigma(xFF_{reco}-xFF_{truth}) (mm)");
	grResolXFF[iv][ipu]->SetLineColor(ipu+1);
	grResolXFF[iv][ipu]->SetMarkerColor(ipu+1);
	grResolXFF[iv][ipu]->SetMarkerStyle(ipu+21);
	//grResolXFF[iv][ipu]->SetMaximum(4);//max_pos);
	//grResolXFF[iv][ipu]->SetMinimum(0.08);//max_pos);
	//grResolXFF[iv][ipu]->GetXaxis()->SetRangeUser(10,200);
	//if (doVsE) grResolXFF[iv][ipu]->GetXaxis()->SetRangeUser(E(10,eta[ieta]),E(200,eta[ieta]));
	grResolXFF[iv][ipu]->Draw(ipu==0?"APE":"PE");
      
	grResolYFF[iv][ipu]->SetTitle(";#eta;#sigma(y_{FF}^{reco}-y_{FF}^{truth}) (mm)");
	if (doVsE) grResolYFF[iv][ipu]->SetTitle(";E (GeV);#sigma(y_{FF}^{reco}-y_{FF}^{truth}) (mm)");
	//grResolYFF[iv][ipu]->SetMaximum(1.);//max_pos);
	//grResolYFF[iv][ipu]->SetMinimum(0.01);//max_pos);
	// grResolYFF[iv][ipu]->GetXaxis()->SetRangeUser(9,110);
	grResolYFF[iv][ipu]->GetYaxis()->SetTitleOffset(0.8);
	//if (doVsE) grResolYFF[iv][ipu]->GetXaxis()->SetRangeUser(E(9,eta[ieta]),E(110,eta[ieta]));
	grResolYFF[iv][ipu]->SetLineColor(ipu+3);
	grResolYFF[iv][ipu]->SetMarkerColor(ipu+3);
	grResolYFF[iv][ipu]->SetMarkerStyle(ipu+23);
	grResolYFF[iv][ipu]->Draw("PE");
      
	mycRP14sig->cd();
	gPad->SetGridy(1);
	//gPad->SetLogy(1);
	grResolX14[iv][ipu]->SetTitle(";#eta;#sigma(x14_{reco}-x14_{truth}) (mm)");
	if (doVsE) grResolX14[iv][ipu]->SetTitle(";E (GeV);#sigma(x14_{reco}-x14_{truth}) (mm)");
	grResolX14[iv][ipu]->SetLineColor(ipu+1);
	grResolX14[iv][ipu]->SetMarkerColor(ipu+1);
	grResolX14[iv][ipu]->SetMarkerStyle(ipu+21);
	//grResolX14[iv][ipu]->SetMaximum(4);//max_pos);
	//grResolX14[iv][ipu]->SetMinimum(0.08);//max_pos);
	//grResolX14[iv][ipu]->GetXaxis()->SetRangeUser(10,200);
	//if (doVsE)  grResolX14[iv][ipu]->GetXaxis()->SetRangeUser(E(10,eta[ieta]),E(200,eta[ieta]));
	grResolX14[iv][ipu]->Draw(ipu==0?"APE":"PE");
      
	grResolY14[iv][ipu]->SetTitle(";#eta;#sigma(y_{14}^{reco}-y_{14}^{truth}) (mm)");
	if (doVsE) grResolY14[iv][ipu]->SetTitle(";E (GeV);#sigma(y_{14}^{reco}-y_{14}^{truth}) (mm)");
	//grResolY14[iv][ipu]->SetMaximum(1.);//max_pos);
	//grResolY14[iv][ipu]->SetMinimum(0.01);//max_pos);
	//grResolY14[iv][ipu]->GetXaxis()->SetRangeUser(9,110);
	grResolY14[iv][ipu]->GetYaxis()->SetTitleOffset(0.8);
	//if (doVsE) grResolY14[iv][ipu]->GetXaxis()->SetRangeUser(E(9,eta[ieta]),E(110,eta[ieta]));
	grResolY14[iv][ipu]->SetLineColor(ipu+3);
	grResolY14[iv][ipu]->SetMarkerColor(ipu+3);
	grResolY14[iv][ipu]->SetMarkerStyle(ipu+23);
	grResolY14[iv][ipu]->Draw("PE");

	mycRAsig->cd();
	gPad->SetGridy(1);
	//gPad->SetLogy(1);
	//grResolTanX[iv][ipu]->SetTitle(";#eta;#sigma(tanA_{reco}-tanA_{truth}) (mrad)");
	//if (doVsE) grResolTanX[iv][ipu]->SetTitle(";E (GeV);#sigma(tanA_{reco}-tanA_{truth}) (mrad)");
	grResolTanX[iv][ipu]->SetTitle(";#eta;#sigma(tanA_{reco}-tanA_{truth}) (mrad)");
	if (doVsE) grResolTanX[iv][ipu]->SetTitle(";E (GeV);#sigma(tanA_{reco}-tanA_{truth}) (mrad)");
	grResolTanX[iv][ipu]->SetLineColor(ipu+1);
	grResolTanX[iv][ipu]->SetMarkerColor(ipu+1);
	grResolTanX[iv][ipu]->SetMarkerStyle(ipu+21);
	//grResolTanX[iv][ipu]->GetXaxis()->SetRangeUser(9,200);
	//if (doVsE) grResolTanX[iv][ipu]->GetXaxis()->SetRangeUser(E(9,eta[ieta]),E(200,eta[ieta]));
	//grResolTanX[iv][ipu]->SetMaximum(10);//max_angle);
	//grResolTanX[iv][ipu]->SetMinimum(1);//max_angle);
	grResolTanX[iv][ipu]->Draw(ipu==0?"APE":"PE");
      
	//grResolTanY[iv][ipu]->SetTitle(";#eta;#sigma(tan(#theta_{y}^{reco})-tan(#theta_{y}^{truth})) (mrad)");
	//if (doVsE) grResolTanY[iv][ipu]->SetTitle(";E (GeV);#sigma((tan(#theta_{y}^{reco})-tan(#theta_{y}^{truth})) (mrad)");
	grResolTanY[iv][ipu]->SetTitle(";#eta;#sigma(#theta_{y}^{reco}-#theta_{y}^{truth}) (mrad)");
	if (doVsE) grResolTanY[iv][ipu]->SetTitle(";E (GeV);#sigma(#theta_{y}^{reco}-#theta_{y}^{truth}) (mrad)");
	//grResolTanY[iv][ipu]->GetXaxis()->SetRangeUser(9,110);
	grResolTanY[iv][ipu]->GetYaxis()->SetTitleOffset(0.8);
	//if (doVsE) grResolTanY[iv][ipu]->GetXaxis()->SetRangeUser(E(9,eta[ieta]),E(110,eta[ieta]));
	//grResolTanY[iv][ipu]->SetMaximum(7);//max_angle);
	//grResolTanY[iv][ipu]->SetMinimum(1);//max_angle);
	grResolTanY[iv][ipu]->SetLineColor(ipu+3);
	grResolTanY[iv][ipu]->SetMarkerColor(ipu+3);
	grResolTanY[iv][ipu]->SetMarkerStyle(ipu+23);
	grResolTanY[iv][ipu]->Draw("PE");

      }//loop on pu

    TLegend *leg1 = new TLegend(0.75,0.75,0.94,0.94);
    leg1->SetFillColor(10);
    leg1->AddEntry(grFit[iv][0],"PU 0","P");
    if (npu>1) leg1->AddEntry(grFit[iv][1],"PU 140","P");
    if (npu==3) leg1->AddEntry(grFit[iv][2],"PU 200","P");

    TLatex lat;
    //lat.SetTextSize(0.04);
    char buf[500];
    //if (eta[ieta]==19){
    mycNice->cd();
    gPad->SetGridy(1);
    for (unsigned ipu(0);ipu<npu;ipu++){
      grResolYFF[iv][ipu]->Draw(ipu==0?"APLE":"PLEsame");
    }
    sprintf(buf,"Single #gamma, p_{T} = %d GeV",pt[0]);
    lat.DrawLatexNDC(0.2,0.85,buf);
    lat.DrawLatexNDC(0.02,0.02,"HGCAL Geant4 Standalone");
    //leg1->Draw("same");
    mycNice->Update();
    if (!doVsE) mycNice->Print(sumDir+"/resolvspT_yFF.C");
    else mycNice->Print(sumDir+"/resolvsE_yFF.C");
    if (!doVsE) mycNice->Print(sumDir+"/resolvspT_yFF.pdf");
    else mycNice->Print(sumDir+"/resolvsE_yFF.pdf");
      
      mycNice->cd();
      gPad->SetGridy(1);
      for (unsigned ipu(0);ipu<npu;ipu++){
	grResolY14[iv][ipu]->Draw(ipu==0?"APLE":"PLEsame");
      }
      lat.DrawLatexNDC(0.02,0.02,"HGCAL Geant4 Standalone");
      lat.DrawLatexNDC(0.2,0.85,buf);
      //leg1->Draw("same");
      mycNice->Update();
      if (!doVsE) mycNice->Print(sumDir+"/resolvspT_y14.C");
      else mycNice->Print(sumDir+"/resolvsE_y14.C");
      if (!doVsE) mycNice->Print(sumDir+"/resolvspT_y14.pdf");
      else mycNice->Print(sumDir+"/resolvsE_y14.pdf");

      mycNice->cd();
      gPad->SetGridy(1);
      for (unsigned ipu(0);ipu<npu;ipu++){
	grResolTanY[iv][ipu]->Draw(ipu==0?"APLE":"PLEsame");
      }
      lat.DrawLatexNDC(0.2,0.85,buf);
      lat.DrawLatexNDC(0.02,0.02,"HGCAL Geant4 Standalone");
      //leg1->Draw("same");
      mycNice->Update();
      if (!doVsE) mycNice->Print(sumDir+"/resolvspT_thetay.C");
      else mycNice->Print(sumDir+"/resolvsE_thetay.C");
      if (!doVsE) mycNice->Print(sumDir+"/resolvspT_tethay.pdf");
      else mycNice->Print(sumDir+"/resolvsE_thetay.pdf");
     



    TLegend *leg2 = new TLegend(0.74,0.66,0.94,0.94);
    leg2->SetFillColor(10);
    leg2->AddEntry(grResidualXFF[iv][0],"PU 0, x","P");
    leg2->AddEntry(grResidualYFF[iv][0],"PU 0, y","P");
    if (npu>1) {
      leg2->AddEntry(grResidualXFF[iv][1],"PU 140, x","P");
      leg2->AddEntry(grResidualYFF[iv][1],"PU 140, y","P");
    }
    if (npu==3) {
      leg2->AddEntry(grResidualXFF[iv][2],"PU 200, x","P");
      leg2->AddEntry(grResidualYFF[iv][2],"PU 200, y","P");
    }

    mycEvt->cd();
    sprintf(buf,"Single #gamma, p_{T}=%d GeV",pt[0]);
    lat.DrawLatexNDC(0.2,0.85,buf);
    lat.DrawLatexNDC(0.02,0.02,"HGCAL Geant4 Standalone");
    leg1->Draw("same");
    
    mycFitQual->cd();
    lat.DrawLatexNDC(0.2,0.85,buf);
    lat.DrawLatexNDC(0.02,0.02,"HGCAL Geant4 Standalone");
    leg1->Draw("same");

    mycRP14->cd();
    lat.DrawLatexNDC(0.2,0.85,buf);
    lat.DrawLatexNDC(0.02,0.02,"HGCAL Geant4 Standalone");
    leg2->Draw("same");

    mycRPFF->cd();
    lat.DrawLatexNDC(0.2,0.85,buf);
    lat.DrawLatexNDC(0.02,0.02,"HGCAL Geant4 Standalone");
    leg2->Draw("same");

    mycRA->cd();
    lat.DrawLatexNDC(0.2,0.85,buf);
    lat.DrawLatexNDC(0.02,0.02,"HGCAL Geant4 Standalone");
    leg2->Draw("same");

    mycRP14sig->cd();
    lat.DrawLatexNDC(0.2,0.85,buf);
    lat.DrawLatexNDC(0.02,0.02,"HGCAL Geant4 Standalone");
    leg2->Draw("same");

    mycRPFFsig->cd();
    lat.DrawLatexNDC(0.2,0.85,buf);
    lat.DrawLatexNDC(0.02,0.02,"HGCAL Geant4 Standalone");
    leg2->Draw("same");

    mycRAsig->cd();
    lat.DrawLatexNDC(0.2,0.85,buf);
    lat.DrawLatexNDC(0.02,0.02,"HGCAL Geant4 Standalone");
    leg2->Draw("same");


  }//loop on pt

  std::ostringstream lPrint;
  
  lPrint.str("");
  lPrint << sumDir << "/SummaryAll_nEvts";
  if (doVsE) lPrint << "_vsE";
  lPrint << ".pdf";
  mycEvt->Update();
  mycEvt->Print(lPrint.str().c_str());
  
  lPrint.str("");
  lPrint << sumDir << "/SummaryAll_fitQuality";
  if (doVsE) lPrint << "_vsE";
  lPrint << ".pdf";
  mycFitQual->Update();
  mycFitQual->Print(lPrint.str().c_str());
  
  lPrint.str("");
  lPrint << sumDir << "/SummaryAll_ResidualPos14";
  if (doVsE) lPrint << "_vsE";
  lPrint << ".pdf";  
  mycRP14->Update();
  mycRP14->Print(lPrint.str().c_str());
  
  lPrint.str("");
  lPrint << sumDir << "/SummaryAll_ResidualPosFF";
  if (doVsE) lPrint << "_vsE";
  lPrint << ".pdf";  
  mycRPFF->Update();
  mycRPFF->Print(lPrint.str().c_str());
  
  lPrint.str("");
  lPrint << sumDir << "/SummaryAll_ResidualAngle";
  if (doVsE) lPrint << "_vsE";
  lPrint << ".pdf";
  mycRA->Update();
  mycRA->Print(lPrint.str().c_str());

  lPrint.str("");
  lPrint << sumDir << "/SummaryAll_ResolPos14";
  if (doVsE) lPrint << "_vsE";
  lPrint << ".pdf";  
  mycRP14sig->Update();
  mycRP14sig->Print(lPrint.str().c_str());
  
  lPrint.str("");
  lPrint << sumDir << "/SummaryAll_ResolPosFF";
  if (doVsE) lPrint << "_vsE";
  lPrint << ".pdf";  
  mycRPFFsig->Update();
  mycRPFFsig->Print(lPrint.str().c_str());
  
  lPrint.str("");
  lPrint << sumDir << "/SummaryAll_ResolAngle";
  if (doVsE) lPrint << "_vsE";
  lPrint << ".pdf";
  mycRAsig->Update();
  mycRAsig->Print(lPrint.str().c_str());

  mycFit->Print(sumDir+"/Positionfits_x.pdf]");
  mycFit->Print(sumDir+"/Positionfits_y.pdf]");

  }//loop on version


  TLatex lat;
  char buf[500];
  TLegend *leg1 = new TLegend(0.75,0.75,0.94,0.94);
  leg1->SetFillColor(10);
  mycNice->cd();
  gPad->SetGridy(1);
  for (unsigned iv(0); iv<nV;++iv){//loop on versions
    for (unsigned ipu(0);ipu<npu;ipu++){
      grResolTanY[iv][ipu]->SetLineColor(iv+1);
      grResolTanY[iv][ipu]->SetMarkerColor(iv+1);
      grResolTanY[iv][ipu]->SetMarkerStyle(iv+21);
      grResolTanY[iv][ipu]->SetMinimum(0);
      grResolTanY[iv][ipu]->SetMaximum(12);
      grResolTanY[iv][ipu]->Draw(iv==0?"APLE":"PLEsame");
    }
  }
  grResolTanY[2][0]->Draw("APLE");
  grResolTanY[1][0]->Draw("PLEsame");
  grResolTanY[0][0]->Draw("PLEsame");
  sprintf(buf,"Single #gamma, p_{T} = %d GeV",pt[0]);
  lat.DrawLatexNDC(0.4,0.85,buf);
  lat.DrawLatexNDC(0.02,0.02,"HGCAL Geant4 Standalone");
  leg1->AddEntry(grResolTanY[0][0],"TP-28","P");
  leg1->AddEntry(grResolTanY[1][0],"TP-24","P");
  leg1->AddEntry(grResolTanY[2][0],"TP-18","P");
  leg1->Draw("same");
  mycNice->Update();
  std::ostringstream lPrint;
  
  lPrint.str("");
  lPrint <<"resolvs";
  if (doVsE) lPrint << "E";
  else lPrint << "eta";
  lPrint << "_thetay_pt" << pt[0] <<".pdf";
  mycNice->Print(lPrint.str().c_str());


  return 0;
    
}//main
