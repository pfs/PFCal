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

int plotResiduals(){

  const unsigned nP = 3;//6;

  unsigned etabin = 19;
  unsigned pt = 60;
  unsigned pu = 0;

  std::ostringstream plotDir[nP];
  plotDir[0] << "/afs/cern.ch/work/a/amagnan/PFCalEEAna/PLOTS/gitV00-02-12/version12/gamma/200um/eta" << etabin << "_et" << pt << "_pu" << pu;
  plotDir[1] << "../PLOTS/gitV00-02-12/version12/gamma/200um/eta" << etabin << "_et" << pt << "_pu" << pu;
  plotDir[2] << "../PLOTS/gitV00-02-12/version12/gamma/200um/eta" << etabin << "_et" << pt << "_pu" << pu;

  std::string suffix[nP] = {
    "",
    "",
    "_new"
    //"_simpleweight_xy",
    //"_simpleweight",
    //"_logweight_xy",
    //"_logweight_all",
    //"_logweight_7_22",
    //"_logweight_linearFrontBack"
  };

  std::string label[nP] = {
    "log weighting, x and y sep",
    "log weighting bug, x and y sep",
    "log weighting 5, x and y sep"
    //"simple weighting, matrix=(x+y)/2.",
    //"simple weighting, x and y sep.",
    //"log weighting, matrix=(x+y)/2.",
    //"log weighting, x and y sep.",
    //"log weighting, fit only l=7-22",
    //"log weighting l=7-22, linear front-back"
  };

  TFile *file[nP];
  SetTdrStyle();
  gStyle->SetOptStat("eMRuo");
  gStyle->SetOptFit(1111);

  TCanvas *mycPx = new TCanvas("mycPx","mycPx",1500,1000);
  mycPx->Divide(3,2);
  TCanvas *mycAx = new TCanvas("mycAx","mycAx",1500,1000);
  mycAx->Divide(3,2);
  TCanvas *mycPy = new TCanvas("mycPy","mycPy",1500,1000);
  mycPy->Divide(3,2);
  TCanvas *mycAy = new TCanvas("mycAy","mycAy",1500,1000);
  mycAy->Divide(3,2);

  TH1F *p_impactX_residual[nP];
  TH1F *p_tanAngleX_residual[nP];
  TH1F *p_impactY_residual[nP];
  TH1F *p_tanAngleY_residual[nP];


  for (unsigned iP(0);iP<nP;++iP){//loop on points
    
    file[iP] = TFile::Open((plotDir[iP].str()+suffix[iP]+".root").c_str());
    if (!file[iP]) {
      std::cout << " -- Error! Cannot open file " << plotDir[iP].str()+suffix[iP] << std::endl;
      return 1;
    }
    file[iP]->cd("PositionFit");

    p_impactX_residual[iP] = (TH1F*)gDirectory->Get("p_impactX14_residual");
    p_tanAngleX_residual[iP] = (TH1F*)gDirectory->Get("p_tanAngleX_residual");
    p_impactY_residual[iP] = (TH1F*)gDirectory->Get("p_impactY14_residual");
    p_tanAngleY_residual[iP] = (TH1F*)gDirectory->Get("p_tanAngleY_residual");
    if (iP==5) {
      p_impactX_residual[iP] = (TH1F*)gDirectory->Get("p_impactX14_residual");
      p_impactY_residual[iP] = (TH1F*)gDirectory->Get("p_impactY14_residual");
    }

    mycPx->cd(iP+1);
    p_impactX_residual[iP]->Draw();
    p_impactX_residual[iP]->Fit("gaus","+","same");
    TLatex lat;
    char buf[500];
    sprintf(buf,"eta=%3.1f, pT=%d GeV, pu=%d",etabin/10.,pt,pu);
    lat.DrawLatexNDC(0.2,0.2,buf);
    lat.DrawLatexNDC(0.2,0.96,label[iP].c_str());

    mycPy->cd(iP+1);
    p_impactY_residual[iP]->Draw();
    p_impactY_residual[iP]->Fit("gaus","+","same");
    lat.DrawLatexNDC(0.2,0.96,label[iP].c_str());
    lat.DrawLatexNDC(0.2,0.2,buf);

    mycAx->cd(iP+1);
    p_tanAngleX_residual[iP]->Draw();
    p_tanAngleX_residual[iP]->Fit("gaus","+","same");
    lat.DrawLatexNDC(0.2,0.96,label[iP].c_str());
    lat.DrawLatexNDC(0.2,0.2,buf);

    mycAy->cd(iP+1);
    p_tanAngleY_residual[iP]->Draw();
    p_tanAngleY_residual[iP]->Fit("gaus","+","same");
    lat.DrawLatexNDC(0.2,0.96,label[iP].c_str());
    lat.DrawLatexNDC(0.2,0.2,buf);

  }//loop on points

  mycPx->Update();
  mycPx->Print((plotDir[1].str()+"_impactXresiduals.pdf").c_str());
  mycPy->Update();
  mycPy->Print((plotDir[1].str()+"_impactYresiduals.pdf").c_str());
  mycAx->Update();
  mycAx->Print((plotDir[1].str()+"_tanAngleXresiduals.pdf").c_str());
  mycAy->Update();
  mycAy->Print((plotDir[1].str()+"_tanAngleYresiduals.pdf").c_str());

  return 0;

}
