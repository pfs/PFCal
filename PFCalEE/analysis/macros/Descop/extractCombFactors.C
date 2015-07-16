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
#include "TProfile.h"

#include "../TDRStyle.h"

int extractCombFactors(){

  unsigned genEnAll[]={3,5,10,30,50};
  const unsigned nGenEnAll=sizeof(genEnAll)/sizeof(unsigned);

  bool useGSC = false;

  const double eta = 2.5;

  const unsigned tpmod = 28;
  std::ostringstream model;
  //if (tpmod>0) model << "_tp" << tpmod;// << "uni";
  if (tpmod>0) model << "_ecal" << tpmod;// << "uni";


  double FHtoEslope = 21.02;
  double FHtoEoffset = 4;
  double BHtoEslope = 8.62;
  double BHtoEoffset = 3.7;
  double ECALslope = 92.83;
  double ECALoffset = 13;
  if (tpmod==2812){
    ECALslope = 93.49;
    ECALoffset = -5;
    FHtoEslope = 22.85;
    FHtoEoffset = 0;//-40;
    BHtoEslope = 4.23;
    BHtoEoffset = 0;//-34.5;
  }
  else if (tpmod==2411){
    ECALslope = 81.64;
    ECALoffset = 2.5;
    FHtoEslope = 22.5;
    FHtoEoffset = 0;//-40;
    BHtoEslope = 5;
    BHtoEoffset = 0;//-34.5;
  }
  else if (tpmod==1809){
    ECALslope = 66.14;
    ECALoffset = 41;
    FHtoEslope = 18.4;
    FHtoEoffset = 0;//-40;
    BHtoEslope = 5;
    BHtoEoffset = 0;//-34.5;
  }
  else if (tpmod==28){
    ECALslope = 92.41;
    ECALoffset = 23;
  }
  else if (tpmod==24){
    ECALslope = 93.35;
    ECALoffset = 14;
    FHtoEslope = 21.59;
    FHtoEoffset = -64;
  }
  else if (tpmod==18){
    ECALslope = 92.26;
    ECALoffset = 6;
    FHtoEslope = 19.08;
    FHtoEoffset = -44;
  }

  TCanvas *mycEtot = new TCanvas("mycEtot","mycEtot",1500,1000);
  mycEtot->Divide(4,2);
  TCanvas *mycEtot1 = new TCanvas("mycEtot1","mycEtot1",1500,1000);
  mycEtot1->Divide(4,2);
  TCanvas *mycSlope = new TCanvas("mycSlope","mycSlope",1500,1000);
  TCanvas *mycEtot2 = new TCanvas("mycEtot2","mycEtot2",1500,1000);
  mycEtot2->Divide(4,2);
  TCanvas *mycEtot3 = new TCanvas("mycEtot3","mycEtot3",1500,1000);
  mycEtot3->Divide(4,2);
  TCanvas *mycSlope2 = new TCanvas("mycSlope2","mycSlope",1500,1000);
  TCanvas *mycSlope3 = new TCanvas("mycSlope3","mycSlope",1500,1000);

  TString plotDir = "/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalDescop/gitV05-02-04/version33/pi-/";

  TFile *inputFile[nGenEnAll];
  TTree *ltree[nGenEnAll];

  TGraphErrors *grSlope = new TGraphErrors();
  TGraphErrors *grSlope2 = new TGraphErrors();
  TGraphErrors *grSlope3 = new TGraphErrors();

  for (unsigned iE(0); iE<nGenEnAll; ++iE){//loop on E
    std::ostringstream linputStr;
    linputStr << plotDir ;
    linputStr << "pion_gc3-50_reso" << model.str() << "_et" << genEnAll[iE] << ".root";
    inputFile[iE] = TFile::Open(linputStr.str().c_str());
    if (!inputFile[iE]) {
      std::cout << " -- Error, input file " << linputStr.str() << " cannot be opened. Skipping..." << std::endl;
      //	    return 1;
      continue;
    }
    else {
      //inputFile[iE]->cd("Energies");
      ltree[iE] = (TTree*)gDirectory->Get("Ereso");
      std::ostringstream lName;

      mycEtot->cd(iE+1);
      lName.str("");
      //lName << "(";
      lName << "(EBHCAL-" << BHtoEoffset << ")/" << BHtoEslope
	//<< ":(EECAL-" << ECALoffset << ")/" << ECALslope
	    << ":(EFHCAL-" << FHtoEoffset << ")/" << FHtoEslope;
      if (useGSC) lName << "*globalC_15";
      //if (genEnAll[iE]<60) lName << "*globalC_5";
      //else if (genEnAll[iE]<300) lName << "*globalC_10";
      //else lName << "*globalC_20";

      std::ostringstream lCut;
      lCut << "EECAL<150";
      //lName << " - " << offset[0][iL][iSR] << ")/" << calib[0][iL][iSR];
      //std::cout << lName.str() << std::endl;
      ltree[iE]->Draw(lName.str().c_str(),lCut.str().c_str(),"colz");

      lName.str("");
      lName << "BHvsFH_" << genEnAll[iE];
      TH2F *hbhfh = (TH2F*)(gPad->GetPrimitive("htemp"))->Clone(lName.str().c_str());
      hbhfh->ProfileX();
      lName << "_pfx";
      TProfile *h2_pfx = (TProfile*)gDirectory->Get(lName.str().c_str());
      h2_pfx->SetMarkerStyle(22);
      h2_pfx->SetMarkerColor(1);
      h2_pfx->Draw("PEsame");
      h2_pfx->Fit("pol1","","same");
      TF1 *fit = (TF1*)h2_pfx->GetFunction("pol1");
      //grSlope->SetPoint(iE,genEnAll[iE],fit->GetParameter(1));
      //grSlope->SetPointError(iE,0,fit->GetParError(1));
      TF1 *fitref = new TF1("fitref","[0]+[1]*x",0,10000);
      fitref->SetParameters(h2_pfx->GetBinContent(1),-1.49);
      if (tpmod==1809) fitref->SetParameters(h2_pfx->GetBinContent(1),-0.99);
      else if (tpmod==2411) fitref->SetParameters(h2_pfx->GetBinContent(1),-1.22);
      fitref->SetLineColor(7);
      fitref->Draw("same");
      grSlope->SetPoint(iE,genEnAll[iE],fabs(fit->GetParameter(1)));
      grSlope->SetPointError(iE,0,fit->GetParError(1));

    }
  }//loop on E
  mycEtot->Update();
  mycEtot->Print(("BHvsFHfits"+model.str()+".pdf").c_str());

  mycSlope->cd();
  gStyle->SetOptFit(1111);
  grSlope->SetTitle(";E_{gen} (GeV);BH/FH");
  grSlope->Draw("APE");
  grSlope->Fit("pol0","","same",40,200);

  TF1 *fit = (TF1*)grSlope->GetFunction("pol0");
  double BHoverFH = fabs(fit->GetParameter(0));


  mycSlope->Update();
  mycSlope->Print(("BHvsFH_slopevsE"+model.str()+".pdf").c_str());

  //return 1;

  double slopeFvsE[nGenEnAll];
  double slopeFvsEerr[nGenEnAll];
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat("e");
  for (unsigned iE(0); iE<nGenEnAll; ++iE){//loop on E

      mycEtot2->cd(iE+1);
      std::ostringstream lName;
      lName.str("");
      //lName << "(";
      lName << "(EBHCAL-" << BHtoEoffset << ")/(" << BHtoEslope << "*" << BHoverFH << ")"
	    << "+(EFHCAL-" << FHtoEoffset << ")/" << FHtoEslope;
      if (useGSC) lName << "*globalC_15";
     //if (genEnAll[iE]<60) lName << "*globalC_5";
      //else if (genEnAll[iE]<300) lName << "*globalC_10";
      //else lName << "*globalC_20";
      lName << ":(EECAL-" << ECALoffset << ")/" << ECALslope;
      //lName << " - " << offset[0][iL][iSR] << ")/" << calib[0][iL][iSR];
      //std::cout << lName.str() << std::endl;
      ltree[iE]->Draw(lName.str().c_str(),"","colz");

      lName.str("");
      lName << "HvsE_" << genEnAll[iE];
      TH2F *hhvse = (TH2F*)(gPad->GetPrimitive("htemp"))->Clone(lName.str().c_str());
      hhvse->Draw("colz");
      hhvse->ProfileX();
      lName << "_pfx";
      TProfile *h2_pfx = (TProfile*)gDirectory->Get(lName.str().c_str());
      h2_pfx->SetMarkerStyle(22);
      h2_pfx->SetMarkerColor(1);
      h2_pfx->Draw("PEsames");
      h2_pfx->Fit("pol1","","same");
      fit = (TF1*)h2_pfx->GetFunction("pol1");

      slopeFvsE[iE] = fabs(fit->GetParameter(1));
      slopeFvsEerr[iE] = fit->GetParError(1);

      grSlope2->SetPoint(iE,genEnAll[iE]*cosh(eta),slopeFvsE[iE]);
      grSlope2->SetPointError(iE,0,slopeFvsEerr[iE]);
      
    }
  mycEtot2->Update();
  mycEtot2->Print(("HvsEfits"+model.str()+".pdf").c_str());

  mycSlope2->cd();
  gStyle->SetOptFit(1111);
  grSlope2->SetTitle(";E_{gen} (GeV);BH+FH/EE");
  grSlope2->Draw("APE");
  TF1 *fitfunc = new TF1("fitfunc","[0]/x+[1]",0,5000);
  grSlope2->Fit("fitfunc","","same");//,0,200);

  //fit = (TF1*)grSlope2->GetFunction("pol0");
  double HoverE = fabs(fitfunc->GetParameter(1));


  mycSlope2->Update();
  mycSlope2->Print(("HvsE_slopevsEgen"+model.str()+".pdf").c_str());



  for (unsigned iE(0); iE<nGenEnAll; ++iE){//loop on E

      mycEtot1->cd(iE+1);
      std::ostringstream lName;
      lName.str("");
      //lName << "(";
      lName << "((EBHCAL-" << BHtoEoffset << ")/(" << BHtoEslope << "*" << BHoverFH << ")"
	    << "+(EFHCAL-" << FHtoEoffset << ")/" << FHtoEslope ;
      if (useGSC) lName << "*globalC_15";
      //if (genEnAll[iE]<60) lName << "*globalC_5";
      //else if (genEnAll[iE]<300) lName << "*globalC_10";
      //else lName << "*globalC_20";
      lName << ")/" << HoverE
	   << "+(EECAL-" << ECALoffset << ")/" << ECALslope;
      //lName << " - " << offset[0][iL][iSR] << ")/" << calib[0][iL][iSR];
      //std::cout << lName.str() << std::endl;
      ltree[iE]->Draw(lName.str().c_str(),"","colz");


      lName.str("");
      lName << "Etot_" << genEnAll[iE];
      TH1F *hetot = (TH1F*)(gPad->GetPrimitive("htemp"))->Clone(lName.str().c_str());
      hetot->Fit("gaus","","same");
      fit = (TF1*)hetot->GetFunction("gaus");

      grSlope3->SetPoint(iE,fit->GetParameter(1),slopeFvsE[iE]);
      grSlope3->SetPointError(iE,fit->GetParError(1),slopeFvsEerr[iE]);
      
  }

  mycEtot1->Update();
  mycEtot1->Print(("Etot"+model.str()+".pdf").c_str());


  mycSlope3->cd();
  gStyle->SetOptFit(1111);
  grSlope3->SetTitle(";E_{tot} (GeV);slope HvsE");
  grSlope3->Draw("APE");
  //grSlope2->SetLineColor(2);
  //grSlope2->SetMarkerColor(2);
  //grSlope2->SetMarkerStyle(22);
  //grSlope2->Draw("PE");
  //TF1 *fitfreq = new TF1("fitfreq","[2]+TMath::Freq((x-[0])/sqrt([1]*x))",0,5000);
  //fitfreq->SetParameters(300,300,-5);
  grSlope3->Fit("fitfunc","","same");//,40,5000);

  fitfunc->SetParameters(30,5.31);
  fitfunc->SetLineColor(4);
  fitfunc->Draw("same");
  //fit = (TF1*)grSlope3->GetFunction("pol0");
  //double HoverE = fabs(fitfunc->GetParameter(1));


  mycSlope3->Update();
  mycSlope3->Print(("HvsE_slopevsEreco"+model.str()+".pdf").c_str());


  return 1;

  const unsigned nGF = 10;

  TCanvas *mycE[nGenEnAll];
  for (unsigned iE(0); iE<nGenEnAll;++iE){
    std::ostringstream lName;
    lName << "mycE" << genEnAll[iE];
    mycE[iE] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1500,1000);
    mycE[iE]->Divide(5,2);
  }


  unsigned list[nGF] = {3,5,10,15,20,25,30,35,40,50};
  for (unsigned iE(0); iE<nGenEnAll; ++iE){//loop on E
    //loop on Cglobal
    for (unsigned ic(0);ic<nGF;++ic){
      mycE[iE]->cd(ic+1);
      std::ostringstream lName;
      lName.str("");
      lName << "((EECAL-" << ECALoffset << ")/" << ECALslope
		    << "+((EFHCAL-" << FHtoEoffset << ")/" << FHtoEslope
		    << "+(EBHCAL-" << BHtoEoffset << ")/(" << BHtoEslope << "*" << BHoverFH << "))/" << HoverE 
	    << "):globalC_" << list[ic];
      ltree[iE]->Draw(lName.str().c_str(),"","colz");
      lName.str("");
      lName << "EvsC_" << genEnAll[iE];
      TH2F *hevsc = (TH2F*)(gPad->GetPrimitive("htemp"))->Clone(lName.str().c_str());
      hevsc->GetXaxis()->SetRangeUser(0.6,1.6);
      hevsc->Draw("colz");
    }
    
    mycE[iE]->Update();
    std::ostringstream lSave;
    lSave << "ErawvsCglobal_" << genEnAll[iE] << model.str() << ".pdf";
    mycE[iE]->Print(lSave.str().c_str());
  }


  return 0;
}//main

