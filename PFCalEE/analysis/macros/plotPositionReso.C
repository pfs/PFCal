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
  double residual_x;
  double residual_y;
  double residual_tanx;
  double residual_tany;
  double residual_x_rms;
  double residual_y_rms;
  double residual_tanx_rms;
  double residual_tany_rms;
};


Result plotOnePoint(unsigned etabin, unsigned pt, unsigned pu){//main

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
  res.residual_x = 0;
  res.residual_y = 0;
  res.residual_tanx = 0;
  res.residual_tany = 0;
  res.residual_x_rms = 0;
  res.residual_y_rms = 0;
  res.residual_tanx_rms = 0;
  res.residual_tany_rms = 0;

  std::ostringstream plotDir;
  plotDir << "../PLOTS/gitV00-02-12/version12/gamma/200um/eta" << etabin << "_et" << pt << "_pu" << pu;

  TFile *file = TFile::Open((plotDir.str()+".root").c_str());
  
  if (!file){
    std::cout << " -- Error, input file " << plotDir.str() << ".root cannot be opened. Skipping..." << std::endl;
    return res;
  }
  else std::cout << " -- input file " << plotDir.str() << ".root successfully opened." << std::endl;

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
    mycE->Print((plotDir.str()+"/AveragePuE"+suffix+".pdf").c_str());
  }
  //return 1;

  TCanvas *mycL = new TCanvas("mycL","mycL",1500,1000);
  TCanvas *mycD = new TCanvas("mycD","mycD",1500,1000);
  TCanvas *mycR = new TCanvas("mycR","mycR",1500,1000);

  TH2D *p_errorMatrix_x = (TH2D*)gDirectory->Get("p_errorMatrix_x");
  TH2D *p_errorMatrix_y = (TH2D*)gDirectory->Get("p_errorMatrix_y");
  TH1F *p_chi2overNDF = (TH1F*)gDirectory->Get("p_chi2overNDF");
  TH1F *p_impactX = (TH1F*)gDirectory->Get("p_impactX");
  TH1F *p_impactX_residual = (TH1F*)gDirectory->Get("p_impactX_residual");
  TH1F *p_impactX_truth = (TH1F*)gDirectory->Get("p_impactX_truth");
  TH1F *p_tanAngleX = (TH1F*)gDirectory->Get("p_tanAngleX");
  TH1F *p_tanAngleX_residual = (TH1F*)gDirectory->Get("p_tanAngleX_residual");
  TH1F *p_tanAngleX_truth = (TH1F*)gDirectory->Get("p_tanAngleX_truth");
  TH1F *p_impactY = (TH1F*)gDirectory->Get("p_impactY");
  TH1F *p_impactY_residual = (TH1F*)gDirectory->Get("p_impactY_residual");
  TH1F *p_impactY_truth = (TH1F*)gDirectory->Get("p_impactY_truth");
  TH1F *p_tanAngleY = (TH1F*)gDirectory->Get("p_tanAngleY");
  TH1F *p_tanAngleY_residual = (TH1F*)gDirectory->Get("p_tanAngleY_residual");
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


  mycL->Divide(3,2);
  mycL->cd(1);
  gPad->SetLogz(1);
  p_errorMatrix_x->SetStats(0);
  p_errorMatrix_x->SetMinimum(0.01);
  p_errorMatrix_x->Draw("colz");

  mycL->cd(4);
  gPad->SetLogz(1);
  p_errorMatrix_y->SetStats(0);
  p_errorMatrix_y->SetMinimum(0.01);
  p_errorMatrix_y->Draw("colz");
  //  gPad->SetLogy(1);
  //gStyle->SetStatW(0.4);
  //gStyle->SetStatH(0.3);
  p_chi2overNDF->GetXaxis()->SetRangeUser(0,20);
  // p_chi2overNDF->Draw();

  res.fitQuality = p_chi2overNDF->GetMean();
  res.fitQuality_rms = p_chi2overNDF->GetRMS();

  mycL->cd(2);
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

  res.nEvts = p_chi2overNDF->GetEntries();
  res.impactX = p_impactX->GetMean();
  res.impactX_rms = p_impactX->GetRMS();

  mycL->cd(3);
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

  mycL->cd(5);
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

  mycL->cd(6);
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
  mycL->Print((plotDir.str()+"/PositionFitSummary"+suffix+".pdf").c_str());

  mycD->Divide(3,2);
  mycD->cd(1);
  gPad->SetLogy(1);
  p_nLayersFit->Draw();
  //p_etavsphi_max->Draw("colz");
  //p_etavsphi_max->SetStats(0);
  mycD->cd(4);
  p_recoZvsLayer->Draw("colz");
  p_recoZvsLayer->SetStats(0);
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
  mycD->Print((plotDir.str()+"/PositionFitDebug"+suffix+".pdf").c_str());

  mycR->Divide(3,2);
  mycR->cd(1);
  p_positionReso->GetXaxis()->SetRangeUser(p_positionReso->GetMean()-5*p_positionReso->GetRMS(),
					    p_positionReso->GetMean()+5*p_positionReso->GetRMS());
  p_positionReso->Draw();
  mycR->cd(4);
  p_angularReso->GetXaxis()->SetRangeUser(0,0.1);
  //p_angularReso->GetMean()-5*p_angularReso->GetRMS(),
  //p_angularReso->GetMean()+5*p_angularReso->GetRMS());
  p_angularReso->Draw();
  mycR->cd(2);
  p_impactX_residual->GetXaxis()->SetRangeUser(p_impactX_residual->GetMean()-5*p_impactX_residual->GetRMS(),
					    p_impactX_residual->GetMean()+5*p_impactX_residual->GetRMS());
  p_impactX_residual->Draw();
  mycR->cd(5);
  p_impactY_residual->GetXaxis()->SetRangeUser(p_impactY_residual->GetMean()-5*p_impactY_residual->GetRMS(),
					    p_impactY_residual->GetMean()+5*p_impactY_residual->GetRMS());
  p_impactY_residual->Draw();
  mycR->cd(3);
  p_tanAngleX_residual->GetXaxis()->SetRangeUser(p_tanAngleX_residual->GetMean()-5*p_tanAngleX_residual->GetRMS(),
					    p_tanAngleX_residual->GetMean()+5*p_tanAngleX_residual->GetRMS());
  p_tanAngleX_residual->Draw();
  mycR->cd(6);
  p_tanAngleY_residual->GetXaxis()->SetRangeUser(p_tanAngleY_residual->GetMean()-5*p_tanAngleY_residual->GetRMS(),
  						 p_tanAngleY_residual->GetMean()+5*p_tanAngleY_residual->GetRMS());
  p_tanAngleY_residual->Draw();

  res.residual_x = p_impactX_residual->GetMean();
  res.residual_x_rms = p_impactX_residual->GetRMS();
  res.residual_y = p_impactY_residual->GetMean();
  res.residual_y_rms = p_impactY_residual->GetRMS();

  res.residual_tanx = p_tanAngleX_residual->GetMean();
  res.residual_tanx_rms = p_tanAngleX_residual->GetRMS();
  res.residual_tany = p_tanAngleY_residual->GetMean();
  res.residual_tany_rms = p_tanAngleY_residual->GetRMS();


  mycR->Update();
  mycR->Print((plotDir.str()+"/PositionFitReso"+suffix+".pdf").c_str());

  return res;

}//main

int plotPositionReso(){//main

  SetTdrStyle();

  std::cout << " Run all ?" << std::endl;
  bool runAll = true;
  std::cin>>runAll;

  const unsigned npu = 1;//2;
  const unsigned npt = 13;

  unsigned pu[npu] = {0};//,140};
  unsigned pt[npt] = {20,30,40,50,60,70,80,90,100,125,150,175,200};
  
  if (!runAll) {
    std::cout << " Please enter pt val:" << std::endl;
    std::cin>>pt[0];
  }

  TCanvas *mycFit = new TCanvas("mycFit","mycFit",1500,1000);
  TCanvas *mycEvt = new TCanvas("mycEvt","mycEvt",1500,1000);
  TCanvas *mycRP = new TCanvas("mycRP","mycRP",1500,1000);
  TCanvas *mycRA = new TCanvas("mycRA","mycRA",1500,1000);
  if (runAll) {
    mycFit->Divide(5,3);
    mycEvt->Divide(5,3);
    mycRP->Divide(5,3);
    mycRA->Divide(5,3);
  }

  TGraphErrors *grEvts[npt][npu];
  TGraphErrors *grFit[npt][npu];
  TGraphErrors *grResidualX[npt][npu];
  TGraphErrors *grResidualY[npt][npu];
  TGraphErrors *grResidualTanX[npt][npu];
  TGraphErrors *grResidualTanY[npt][npu];

  const unsigned neta = 1;//7;
  const unsigned eta[neta] = {17};//,19,21,23,25,27,29};
  double etaval[npu][neta];

  
  for (unsigned ipt=0; ipt<(runAll?npt:1);++ipt){//loop on pt values
    for (unsigned ipu(0); ipu<npu;++ipu){//loop on pu values
      
      
      std::vector<Result> resVec;

      unsigned newp = 0;
      for (unsigned ie(0);ie<neta;++ie){
	Result lres = plotOnePoint(eta[ie],pt[ipt],pu[ipu]);
	if (lres.nEvts > 0){
	  etaval[ipu][newp] = eta[ie]/10.+0.02*(2.*ipu-1.);
	  //std::cout << pu[ipu] << " " << eta[ie] << " " << etaval[ipu][newp] << std::endl;
	  newp++;
	  resVec.push_back(lres);
	}
	if (!runAll) {
	  /*std::cout << " Proceed to next point ?" << std::endl;
	  bool lNext = 0;
	  std::cin>>lNext;
	
	  if (lNext) continue;
	  else break;*/
	}
      }
      //continue;

      const unsigned nP = resVec.size();
      double etaerr[nP];
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
      double residual_x[nP];
      double residual_y[nP];
      double residual_tanx[nP];
      double residual_tany[nP];
      double residual_x_rms[nP];
      double residual_y_rms[nP];
      double residual_tanx_rms[nP];
      double residual_tany_rms[nP];
      
      double max_angle = 0;
      double max_pos = 0;
      for (unsigned iP(0); iP<nP;++iP){
	etaerr[iP] = 0;
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
	residual_x[iP] = resVec[iP].residual_x;
	residual_y[iP] = resVec[iP].residual_y;
	residual_tanx[iP] = resVec[iP].residual_tanx*1000;
	residual_tany[iP] = resVec[iP].residual_tany*1000;
	residual_x_rms[iP] = resVec[iP].residual_x_rms;
	residual_y_rms[iP] = resVec[iP].residual_y_rms;
	residual_tanx_rms[iP] = resVec[iP].residual_tanx_rms*1000;
	residual_tany_rms[iP] = resVec[iP].residual_tany_rms*1000;
	if ( (fabs(residual_tanx[iP])+residual_tanx_rms[iP])> max_angle) max_angle = fabs(residual_tanx[iP])+residual_tanx_rms[iP];
	if ( (fabs(residual_tany[iP])+residual_tany_rms[iP])> max_angle) max_angle = fabs(residual_tany[iP])+residual_tany_rms[iP];
	if ( (fabs(residual_x[iP])+residual_x_rms[iP])> max_pos) max_pos = fabs(residual_x[iP])+residual_x_rms[iP];
	if ( (fabs(residual_y[iP])+residual_y_rms[iP])> max_pos) max_pos = fabs(residual_y[iP])+residual_y_rms[iP];
      }
      
      grEvts[ipt][ipu] = new TGraphErrors(nP,etaval[ipu],nEvts,etaerr,etaerr);
      grFit[ipt][ipu] = new TGraphErrors(nP,etaval[ipu],fitQuality,etaerr,fitQuality_rms);

      grResidualX[ipt][ipu] = new TGraphErrors(nP,etaval[ipu],residual_x,etaerr,residual_x_rms);
      grResidualY[ipt][ipu] = new TGraphErrors(nP,etaval[ipu],residual_y,etaerr,residual_y_rms);
      grResidualTanX[ipt][ipu] = new TGraphErrors(nP,etaval[ipu],residual_tanx,etaerr,residual_tanx_rms);
      grResidualTanY[ipt][ipu] = new TGraphErrors(nP,etaval[ipu],residual_tany,etaerr,residual_tany_rms);

      if (runAll) mycEvt->cd(ipt+1);
      else mycEvt->cd();
      gPad->SetGridy(1);
      grEvts[ipt][ipu]->SetTitle(";#eta;N_{events}");
      grEvts[ipt][ipu]->SetLineColor(ipu+1);
      grEvts[ipt][ipu]->SetMarkerColor(ipu+1);
      grEvts[ipt][ipu]->SetMarkerStyle(ipu+21);
      grEvts[ipt][ipu]->SetMinimum(0);
      grEvts[ipt][ipu]->SetMaximum(1200);
      grEvts[ipt][ipu]->Draw(ipu==0?"APE":"PE");
      
      if (runAll) mycFit->cd(ipt+1);
      else mycFit->cd();
      gPad->SetGridy(1);
      grFit[ipt][ipu]->SetTitle(";#eta;#chi^{2}/NDF");
      grFit[ipt][ipu]->SetLineColor(ipu+1);
      grFit[ipt][ipu]->SetMarkerColor(ipu+1);
      grFit[ipt][ipu]->SetMarkerStyle(ipu+21);
      grFit[ipt][ipu]->SetMinimum(0);
      grFit[ipt][ipu]->SetMaximum(10);
      grFit[ipt][ipu]->Draw(ipu==0?"APE":"PE");
      
      if (runAll) mycRP->cd(ipt+1);
      else mycRP->cd();
      grResidualX[ipt][ipu]->SetTitle(";#eta;x_{reco}-x_{truth} (mm)");
      grResidualX[ipt][ipu]->SetLineColor(ipu+1);
      grResidualX[ipt][ipu]->SetMarkerColor(ipu+1);
      grResidualX[ipt][ipu]->SetMarkerStyle(ipu+21);
      grResidualX[ipt][ipu]->SetMaximum(4);//max_pos);
      grResidualX[ipt][ipu]->SetMinimum(-4);//max_pos);
      grResidualX[ipt][ipu]->Draw(ipu==0?"APE":"PE");
      
      grResidualY[ipt][ipu]->SetLineColor(ipu+3);
      grResidualY[ipt][ipu]->SetMarkerColor(ipu+3);
      grResidualY[ipt][ipu]->SetMarkerStyle(ipu+23);
      grResidualY[ipt][ipu]->Draw("PE");
      
      if (runAll) mycRA->cd(ipt+1);
      else mycRA->cd();
      grResidualTanX[ipt][ipu]->SetTitle(";#eta;tanA_{reco}-tanA_{truth} (mrad)");
      grResidualTanX[ipt][ipu]->SetLineColor(ipu+1);
      grResidualTanX[ipt][ipu]->SetMarkerColor(ipu+1);
      grResidualTanX[ipt][ipu]->SetMarkerStyle(ipu+21);
      grResidualTanX[ipt][ipu]->SetMaximum(40);//max_angle);
      grResidualTanX[ipt][ipu]->SetMinimum(-40);//max_angle);
      grResidualTanX[ipt][ipu]->Draw(ipu==0?"APE":"PE");
      
      grResidualTanY[ipt][ipu]->SetLineColor(ipu+3);
      grResidualTanY[ipt][ipu]->SetMarkerColor(ipu+3);
      grResidualTanY[ipt][ipu]->SetMarkerStyle(ipu+23);
      grResidualTanY[ipt][ipu]->Draw("PE");
      
    }//loop on pu

    if (npu>1){
      TLegend *leg1 = new TLegend(0.64,0.8,0.94,0.94);
      leg1->SetFillColor(10);
      leg1->AddEntry(grFit[ipt][0],"PU 0","P");
      leg1->AddEntry(grFit[ipt][1],"PU 140","P");
      
      TLegend *leg2 = new TLegend(0.64,0.66,0.94,0.94);
      leg2->SetFillColor(10);
      leg2->AddEntry(grResidualX[ipt][0],"PU 0, x","P");
      leg2->AddEntry(grResidualY[ipt][0],"PU 0, y","P");
      leg2->AddEntry(grResidualX[ipt][1],"PU 140, x","P");
      leg2->AddEntry(grResidualY[ipt][1],"PU 140, y","P");
      
      if (runAll) mycEvt->cd(ipt+1);
      else mycEvt->cd();
      TLatex lat;
      //lat.SetTextSize(0.04);
      char buf[500];
      sprintf(buf,"E_{T}^{#gamma} = %d GeV",pt[ipt]);
      lat.DrawLatexNDC(0.2,0.85,buf);
      lat.DrawLatexNDC(0.02,0.02,"HGCAL Geant4 Standalone");
      leg1->Draw("same");
      
      if (runAll) mycFit->cd(ipt+1);
      else mycFit->cd();
      lat.DrawLatexNDC(0.2,0.85,buf);
      lat.DrawLatexNDC(0.02,0.02,"HGCAL Geant4 Standalone");
      leg1->Draw("same");
      
      if (runAll) mycRP->cd(ipt+1);
      else mycRP->cd();
      lat.DrawLatexNDC(0.2,0.85,buf);
      lat.DrawLatexNDC(0.02,0.02,"HGCAL Geant4 Standalone");
      leg2->Draw("same");
      
      if (runAll) mycRA->cd(ipt+1);
      else mycRA->cd();
      lat.DrawLatexNDC(0.2,0.85,buf);
      lat.DrawLatexNDC(0.02,0.02,"HGCAL Geant4 Standalone");
      leg2->Draw("same");
    }

  }//loop on pt
  
  
  std::ostringstream lPrint;
  
  lPrint.str("");
  lPrint << "PLOTS/Summary_nEvts";
  if (!runAll) lPrint << "_et" << pt[0];
  lPrint << ".pdf";
  mycEvt->Update();
  mycEvt->Print(lPrint.str().c_str());
    
  lPrint.str("");
  lPrint << "PLOTS/Summary_fitQuality";
  if (!runAll) lPrint << "_et" << pt[0];
  lPrint << ".pdf";
  mycFit->Update();
  mycFit->Print(lPrint.str().c_str());
  
  lPrint.str("");
  lPrint << "PLOTS/Summary_ResidualPos";
  if (!runAll) lPrint << "_et" << pt[0];
  lPrint << ".pdf";  
  mycRP->Update();
  mycRP->Print(lPrint.str().c_str());
  
  lPrint.str("");
  lPrint << "PLOTS/Summary_ResidualAngle";
  if (!runAll) lPrint << "_et" << pt[0];
  lPrint << ".pdf";
  mycRA->Update();
  mycRA->Print(lPrint.str().c_str());

  return 0;

}//main
