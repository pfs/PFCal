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

int plotResovsRemove(){

  const unsigned nRemove = 13;
  //unsigned lToRemove[nRemove] = {27,24,21,8,19,3,6,1,10,14,17,12};

  int colorSel[nRemove] = {kBlack,kBlue+2,kBlue,kBlue-2,kCyan,
			   kGreen+2,kGreen,kGreen-2,kYellow+2,kOrange+2,
			   kPink,kRed-2,kRed};
  int color[nRemove] = {kGray,kBlue+1,kBlue-1,kBlue-3,kCyan-1,
			kGreen+1,kGreen-1,kGreen-3,kYellow+1,kOrange+1,
			kPink-1,kRed-3,kRed-1};

  TF1 *fit[nRemove];
  TF1 *fitRef[nRemove];

  const unsigned SR = 4;
  const unsigned nL = 31;
  TGraphErrors *grRes[nRemove][nL];
  TGraphErrors *grCal[nRemove][nL];
  TGraphErrors *grDelta[nRemove][nL];
  TGraphErrors* grOffset[nRemove];
  TGraphErrors* grSlope[nRemove];

  double id[nRemove];
  double did[nRemove];
  double s[3][nRemove];
  double c[3][nRemove];
  double ds[3][nRemove];
  double dc[3][nRemove];

  const unsigned nCanvas = 2;//nRemove;  
  TCanvas *mycC[nCanvas];
  TCanvas *mycR[nCanvas];
  TCanvas *mycS[nCanvas];
  TCanvas *mycF[nCanvas];
  for (unsigned iC(0);iC<nCanvas;++iC){
    std::ostringstream lName;
    lName << "mycC" << iC;
    if (iC<1) mycC[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1500,1000);
    if (iC<1) mycC[iC]->Divide(2,2);
    lName.str("");
    lName << "mycR" << iC;
    if (iC<1) mycR[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1500,1000);
    lName.str("");
    lName << "mycF" << iC;
    mycF[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1500,1000);
    lName.str("");
    lName << "mycS" << iC;
    if (iC<1) mycS[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1500,1000);
    if (iC<1) mycS[iC]->Divide(2,2);
  }

  TLegend *leg = new TLegend(0.85,0.42,0.99,0.99);
  leg->SetFillColor(10);
  TLegend *legF = new TLegend(0.85,0.42,0.99,0.99);
  legF->SetFillColor(10);
  TLegend *legFR = new TLegend(0.85,0.62,0.99,0.99);
  legFR->SetFillColor(10);

  TFile *fin[nRemove];
  bool first = true;
  for (unsigned r(0); r<nRemove;++r){//loop on nremove
    id[r] = r;
    did[r] = 0;
    std::ostringstream lname;
    lname << "./PLOTS/CalibReso_vsE_Remove";
    if (r<12) lname << r ;
    else lname << r-1 ;
    lname <<"_nofit.root";
    fin[r] = TFile::Open(lname.str().c_str());
    for (unsigned l(0); l<nL;++l){//loop on layers
      bool selected = (r<12 && l==30) || (r==12 && l==12);
      lname.str("");
      lname << "SR" << SR << "_remove" << l;
      fin[r]->cd(lname.str().c_str());
      lname.str("");
      lname << "calibRecoFit" << SR << "lay" << l << "pu0";
      grCal[r][l] = (TGraphErrors*)gDirectory->Get(lname.str().c_str());
      mycC[0]->cd(1);
      grCal[r][l]->SetLineColor(selected ? colorSel[r]: color[r]);
      grCal[r][l]->SetMarkerColor(selected ? colorSel[r]: color[r]);
      grCal[r][l]->SetMarkerStyle(r+20);
      grCal[r][l]->GetXaxis()->SetTitle("E_{gen} (GeV)");
      grCal[r][l]->GetXaxis()->SetLabelSize(0.05);
      grCal[r][l]->GetXaxis()->SetTitleSize(0.05);
      grCal[r][l]->GetXaxis()->SetTitleOffset(1);
      grCal[r][l]->GetYaxis()->SetLabelSize(0.035);
      grCal[r][l]->GetYaxis()->SetTitleSize(0.05);
      grCal[r][l]->GetYaxis()->SetTitleOffset(1);
      if (selected) grCal[r][l]->SetLineWidth(3);
      if (selected) grCal[r][l]->Draw(first?"APEL":"PEL");
      lname.str("");
      lname << "calibRecoDelta" << SR << "lay" << l << "pu0";
      grDelta[r][l] = (TGraphErrors*)gDirectory->Get(lname.str().c_str());
      mycC[0]->cd(3);
      grDelta[r][l]->SetLineColor(selected ? colorSel[r]: color[r]);
      grDelta[r][l]->SetMarkerColor(selected ? colorSel[r]: color[r]);
      grDelta[r][l]->SetMarkerStyle(r+20);
      grDelta[r][l]->GetYaxis()->SetRangeUser(-0.015,0.015);
      grDelta[r][l]->GetXaxis()->SetLabelSize(0.05);
      grDelta[r][l]->GetXaxis()->SetTitleSize(0.05);
      grDelta[r][l]->GetXaxis()->SetTitleOffset(1);
      grDelta[r][l]->GetYaxis()->SetLabelSize(0.05);
      grDelta[r][l]->GetYaxis()->SetTitleSize(0.05);
      grDelta[r][l]->GetYaxis()->SetTitleOffset(1);
      if (selected) grDelta[r][l]->SetLineWidth(3);
      if (selected) grDelta[r][l]->Draw(first?"APEL":"PEL");
      lname.str("");
      lname << "resoRecoFit" << SR << "lay" << l << "pu1";
      grRes[r][l] = (TGraphErrors*)gDirectory->Get(lname.str().c_str());
      mycR[0]->cd();
      grRes[r][l]->SetLineColor(selected ? colorSel[r]: color[r]);
      grRes[r][l]->SetMarkerColor(selected ? colorSel[r]: color[r]);
      grRes[r][l]->SetMarkerStyle(r+20);
      //if (selected) grRes[r][l]->SetLineWidth(3);
      if (selected) {
	grRes[r][l]->Draw(first?"APEL":"PEL");
	lname.str("");
	lname << "-" << r << " lay";
	leg->AddEntry(grRes[r][l],lname.str().c_str(),"PL");
      }
      if (selected) {
	lname.str("");
	lname << "reso_" << r;
	fit[r] = new TF1(lname.str().c_str(),"sqrt([0]*[0]/x+[1]*[1]+[2]*[2]/(x*x))",0,1000);
	fit[r]->SetParameters(0.212,0.009,0.04);
	fit[r]->FixParameter(2,0.04);
	grRes[r][l]->Fit(fit[r],"BIME0Q");
	std::cout << "Fit" << r << " " <<  fit[r]->GetParameter(0) << " " << fit[r]->GetParameter(1) << std::endl;

	fit[r]->SetLineColor(colorSel[r]);
	fit[r]->SetLineWidth(2);
	//fit[r]->Draw("same");
	s[0][r] = fit[r]->GetParameter(0);
	ds[0][r] = fit[r]->GetParError(0);
	c[0][r] = fit[r]->GetParameter(1);
	dc[0][r] = fit[r]->GetParError(1);
	s[1][r] = fit[r]->Eval(50);
	ds[1][r] = 0;
 	s[2][r] = fit[r]->Eval(250);
	ds[2][r] = 0;
     }
      
      if (selected) first = false;
    }//loop on layers


    mycR[0]->cd();
    lname.str("");
    lname << "resoRef_" << r;
    fitRef[r] = new TF1(lname.str().c_str(),"[3]*sqrt([0]*[0]/x+[1]*[1]+[2]*[2]/(x*x))",0,1000);
    fitRef[r]->SetParameters(0.212,0.009,0,1);
    fitRef[r]->SetLineColor(1);
    fitRef[r]->SetLineStyle(7);
    fitRef[r]->SetParameter(3,sqrt(30./(30-r)));
    fitRef[r]->SetLineColor(colorSel[r]);
    //fitRef[r]->Draw("same");

    mycC[0]->cd(2);
    lname.str("");
    lname << "SR" << SR ;
    fin[r]->cd(lname.str().c_str());
    lname.str("");
    lname << "grOffset_0_SR" << SR ;
    grOffset[r] = (TGraphErrors*)gDirectory->Get(lname.str().c_str());
    grOffset[r]->SetLineColor(color[r]);
    grOffset[r]->SetMarkerColor(color[r]);
    grOffset[r]->SetMarkerStyle(r+20);
    grOffset[r]->SetMinimum(-30);
    grOffset[r]->SetMaximum(50);
    grOffset[r]->Draw(r==0?"APEL":"PEL");
    mycC[0]->cd(4);
    lname.str("");
    lname << "grSlope_0_SR" << SR ;
    grSlope[r] = (TGraphErrors*)gDirectory->Get(lname.str().c_str());
    grSlope[r]->SetLineColor(color[r]);
    grSlope[r]->SetMarkerColor(color[r]);
    grSlope[r]->SetMarkerStyle(r+20);
    grSlope[r]->SetMinimum(80);
    grSlope[r]->SetMaximum(92);
    grSlope[r]->Draw(r==0?"APEL":"PEL");


  }//loop on nremove

  mycC[0]->Update();
  std::ostringstream save;
  save << "PLOTS/RemoveScenarios_Calibration.pdf";
  mycC[0]->Print(save.str().c_str());

  mycR[0]->cd();

  TFile *fref[3];
  fref[0] = TFile::Open("PLOTS/CalibReso_eta16_vsE_IC3_TP28.root");
  fref[1] = TFile::Open("PLOTS/CalibReso_eta16_vsE_IC3_TP24.root");
  fref[2] = TFile::Open("PLOTS/CalibReso_eta16_vsE_IC3_TP18.root");
  TGraphErrors *grRef[3];
  TF1 *fittp[3];
  std::ostringstream lname;
  for (unsigned ir(0);ir<3;++ir){
    lname.str("");
    lname << "SR" << SR;
    fref[ir]->cd(lname.str().c_str());
    lname.str("");
    lname << "resoRecoFit" << SR << "eta16pu1";
    grRef[ir] = (TGraphErrors*)gDirectory->Get(lname.str().c_str());
    grRef[ir]->SetMarkerStyle(2+ir);
    grRef[ir]->SetMarkerColor(6);
    grRef[ir]->SetLineColor(6);
    grRef[ir]->SetLineStyle(3+ir);
    //grRef[ir]->Draw("PE");
    lname.str("");
    lname << "fitTP_" << ir;
    fittp[ir] = new TF1(lname.str().c_str(),"sqrt([0]*[0]/x+[1]*[1]+[2]*[2]/(x*x))",0,1000);
    fittp[ir]->SetParameters(0.206,0.008,0.04);
    fittp[ir]->FixParameter(2,0.04);
    fittp[ir]->FixParameter(1,0.009);
    fittp[ir]->SetLineColor(6);
    fittp[ir]->SetLineStyle(3+ir);
    grRef[ir]->Fit(fittp[ir],"BIME0");
    fittp[ir]->Draw("same");
    std::cout << ir << " " <<  fittp[ir]->GetParameter(0) << " " << fittp[ir]->GetParameter(1) <<  " " << fittp[ir]->GetParameter(2) << std::endl;
    std::cout << ir << " " <<  fittp[ir]->GetParError(0) << " " << fittp[ir]->GetParError(1) <<  " " << fittp[ir]->GetParError(2) << std::endl;
  }

  leg->AddEntry(fittp[0],"TP-28","L");
  leg->AddEntry(fittp[1],"TP-24","L");
  leg->AddEntry(fittp[2],"TP-18","L");
  leg->Draw("same");

  mycR[0]->Update();
  save.str("");
  save << "PLOTS/RemoveScenarios_Resolution.pdf";
  mycR[0]->Print(save.str().c_str());

  grRes[0][30]->GetXaxis()->SetRangeUser(0,50);
  save.str("");
  save << "PLOTS/RemoveScenarios_Resolution_zoom.pdf";
  mycR[0]->Print(save.str().c_str());

  mycF[0]->cd();
  gPad->SetLogx(1);
  grRes[0][30]->GetXaxis()->SetRangeUser(0,2000);
  grRes[0][30]->Draw("APE");
  fit[0]->Draw("same");
  grRes[2][30]->Draw("PE");
  fit[2]->Draw("same");
  grRes[5][30]->Draw("PE");
  fit[5]->Draw("same");
  grRes[12][12]->Draw("PE");
  fit[12]->Draw("same");
  for (unsigned ir(0);ir<3;++ir){
    grRef[ir]->Draw("PE");
    fittp[ir]->Draw("same");
  }
  legF->AddEntry(grRes[0][30],"CMSSW-30","PL");
  legF->AddEntry(grRes[2][30],"CMSSW-28","PL");
  legF->AddEntry(grRes[5][30],"CMSSW-25","PL");
  legF->AddEntry(grRes[12][12],"CMSSW-18","PL");
  legF->AddEntry(grRef[0],"TP-28","PL");
  legF->AddEntry(grRef[1],"TP-24","PL");
  legF->AddEntry(grRef[2],"TP-18","PL");
  legF->Draw("same");

  mycF[0]->Update();
  save.str("");
  save << "PLOTS/RemoveScenarios_Resolution_final.pdf";
  mycF[0]->Print(save.str().c_str());


  mycF[1]->cd();
  gPad->SetLogx(1);
  gPad->SetGridy(1);
  grRes[0][30]->GetXaxis()->SetRangeUser(0,2000);
  TGraphErrors* grRatios[4];
  grRatios[0] = (TGraphErrors*)grRes[0][30]->Clone("grRatio_30_28");
  grRatios[1] = (TGraphErrors*)grRes[2][30]->Clone("grRatio_28_28");
  grRatios[2] = (TGraphErrors*)grRes[5][30]->Clone("grRatio_25_28");
  grRatios[3] = (TGraphErrors*)grRes[12][12]->Clone("grRatio_18_28");
  grRatios[0]->SetMarkerColor(kGray);
  grRatios[0]->SetLineColor(kGray);
  grRatios[0]->SetMarkerStyle(20);
  grRatios[1]->SetMarkerColor(1);
  grRatios[1]->SetLineColor(1);
  grRatios[1]->SetMarkerStyle(20);
  grRatios[2]->SetMarkerColor(2);
  grRatios[2]->SetLineColor(2);
  grRatios[2]->SetMarkerStyle(21);
  grRatios[3]->SetMarkerColor(4);
  grRatios[3]->SetLineColor(4);
  grRatios[3]->SetMarkerStyle(22);
  grRatios[0]->SetTitle(";E^{GEN}_{#gamma};(#sigma/#mu)/(#sigma/#mu)_{28}");
  for (unsigned iP(0); iP<grRes[2][30]->GetN();++iP){
    double xn,yn,xd,yd;
    double yne = grRes[0][30]->GetErrorY(iP);
    double yde = grRes[2][30]->GetErrorY(iP);
    grRes[0][30]->GetPoint(iP,xn,yn);
    grRes[2][30]->GetPoint(iP,xd,yd);
    double err = yn/yd*sqrt(pow(yne/yn,2)+pow(yde/yd,2));
    grRatios[0]->SetPoint(iP,xn,yn/yd);
    grRatios[0]->SetPointError(iP,0,err);

    std::cout << iP << " 30 " << yn/yd << " " << err << std::endl;

    grRes[2][30]->GetPoint(iP,xn,yn);
    yne = grRes[2][30]->GetErrorY(iP);
    err = yn/yd*sqrt(pow(yne/yn,2)+pow(yde/yd,2));
    grRatios[1]->SetPoint(iP,xn,yn/yd);
    grRatios[1]->SetPointError(iP,0,err);
    std::cout << iP << " 28 " << yn/yd << " " << err << std::endl;

    grRes[5][30]->GetPoint(iP,xn,yn);
    yne = grRes[5][30]->GetErrorY(iP);
    err = yn/yd*sqrt(pow(yne/yn,2)+pow(yde/yd,2));
    grRatios[2]->SetPoint(iP,xn,yn/yd);
    grRatios[2]->SetPointError(iP,0,err);
    std::cout << iP << " 25 " << yn/yd << " " << err << std::endl;

    grRes[12][12]->GetPoint(iP,xn,yn);
    yne = grRes[12][12]->GetErrorY(iP);
    err = yn/yd*sqrt(pow(yne/yn,2)+pow(yde/yd,2));
    grRatios[3]->SetPoint(iP,xn,yn/yd);
    grRatios[3]->SetPointError(iP,0,err);
    std::cout << iP << " 18 " << yn/yd << " " << err << std::endl;
  }

  grRatios[0]->GetYaxis()->SetRangeUser(0.75,1.5);
  grRatios[0]->Draw("APE");
  grRatios[1]->Draw("PE");
  grRatios[2]->Draw("PE");
  grRatios[3]->Draw("PE");
  legFR->AddEntry(grRatios[0],"CMSSW-30/28","P");
  legFR->AddEntry(grRatios[1],"CMSSW-28/28","P");
  legFR->AddEntry(grRatios[2],"CMSSW-25/28","P");
  legFR->AddEntry(grRatios[3],"CMSSW-18/28","P");
  legFR->Draw("same");

  mycF[1]->Update();
  save.str("");
  save << "PLOTS/RemoveScenarios_ResolutionRatios_final.pdf";
  mycF[1]->Print(save.str().c_str());



  TGraphErrors *grSampl = new TGraphErrors(nRemove,id,s[0],did,ds[0]);
  TGraphErrors *grConst = new TGraphErrors(nRemove,id,c[0],did,dc[0]);
  TGraphErrors *grRes50 = new TGraphErrors(nRemove,id,s[1],did,ds[1]);
  TGraphErrors *grRes250 = new TGraphErrors(nRemove,id,s[2],did,ds[2]);

  mycS[0]->cd(1);
  gPad->SetGridy(1);
  grSampl->SetTitle(";# removed;sampling (GeV^{0.5})");
  grSampl->SetMinimum(0.1);
  grSampl->SetMaximum(0.3);
  grSampl->Draw("APL");

  TF1 *ref = new TF1("ref","[0]*sqrt(30./(30.-x))",0,12);
  ref->SetParameter(0,s[0][0]);
  ref->SetLineColor(6);
  ref->Draw("same");

  mycS[0]->cd(2);
  gPad->SetGridy(1);
  grConst->SetTitle(";# removed;Constant");
  grConst->SetMinimum(0);
  grConst->SetMaximum(0.03);
  grConst->Draw("APL");

  mycS[0]->cd(3);
  gPad->SetGridy(1);
  grRes50->SetTitle(";# removed;#sigma/E @ 50 GeV");
  grRes50->SetMinimum(0);
  grRes50->SetMaximum(0.06);
  grRes50->Draw("APL");

  TF1 *ref50 = new TF1("ref","[0]*sqrt(30./(30.-x))",0,12);
  ref50->SetParameter(0,s[1][0]);
  ref50->SetLineColor(6);
  ref50->Draw("same");

  mycS[0]->cd(4);
  gPad->SetGridy(1);
  grRes250->SetTitle(";# removed;#sigma/E @ 250 GeV");
  grRes250->SetMinimum(0);
  grRes250->SetMaximum(0.04);
  grRes250->Draw("APL");

  TF1 *ref250 = new TF1("ref","[0]*sqrt(30./(30.-x))",0,12);
  ref250->SetParameter(0,s[2][0]);
  ref250->SetLineColor(6);
  ref250->Draw("same");


  mycS[0]->Update();
  save.str("");
  save << "PLOTS/RemoveScenarios_compa.pdf";
  mycS[0]->Print(save.str().c_str());


  return 0;

}//main
