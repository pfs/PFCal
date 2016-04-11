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
#include "TPaveStats.h"

#include "effSigmaMacro.C"

void drawLayerWithGap(TLatex & lat, bool random, double ypos=0)
{
  lat.SetTextColor(9);
  lat.SetTextSize(0.03);
  lat.SetTextAngle(90);
  char buf[500];

  for (unsigned il(0);il<28;++il){
    if (random) sprintf(buf,"Layer %d",il);
    else sprintf(buf,"Layer %d-%d",il,il+1);
    double offset = random ? 10.*((7*il)%31) : il/2*30;
    if (-230+offset < -170) continue;
    if (random || (!random && il%2==0)) lat.DrawLatex(-230.+offset,ypos,buf);
  }
  
  for (unsigned il(0);il<28;++il){
    if (random) sprintf(buf,"Layer %d",il);
    else sprintf(buf,"Layer %d-%d",il,il+1);
    double offset = random ? 10.*((7*il)%31) : il/2*30;
    if (230+offset>340) continue;
    if (random || (!random && il%2==0)) lat.DrawLatex(230.+offset,ypos,buf);
  }
  lat.SetTextAngle(0);
  lat.SetTextColor(1);
  lat.SetTextSize(0.05);
}

unsigned getGapLayer(const double & vtxx){

  for (unsigned l(0);l<28;++l){
    double shift = -160+10*((7*l)%31);
    if (vtxx>145) shift = 150+10*((7*l)%31);
    if (fabs(vtxx-shift)<=5) return l;
  }

  return 28;
};

int plotSigmaEffvsVtxPos(){//main


  const unsigned nF = 2;
  TFile *fcalib[nF];
  fcalib[0] = TFile::Open("/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalCracks/gitV06-03-04/version100/model3/gamma/eta20_et60_pu0_IC3_Si2.root");
  fcalib[1] = TFile::Open("/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalCracks/gitV06-03-04/version100/model3/gamma/phi_0.440pi/eta20_et60_pu0_IC3_Si2.root");

  const double eta = 2.0;
  const double pT = 60;
  const double Egen = pT*cosh(eta);

  const std::string version = "V06d-04-06";

  std::string label[nF] = {"phi90","phi79"};

  bool linedup = (version=="V06-03-04");
  bool random = (version=="V06a-03-04");
  if (version=="V06a-03-04") {
    label[0] = "V06aphi90";
    label[1] = "V06aphi79";
  }
  else if (version=="V06b-03-05") {
    label[0] = "V06bphi90";
    label[1] = "V06bphi79";
  }
  else if (version=="V06c-04-06") {
    label[0] = "V06cphi90";
    label[1] = "V06cphi79";
  }
  else if (version=="V06d-04-06") {
    label[0] = "V06dphi90";
    label[1] = "V06dphi79";
  }


  TFile *fin[nF];
  fin[0] = TFile::Open(("/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalCracks/git"+version+"/version100/model4/gamma/eta20_et60_pu0_IC3_Si2.root").c_str());
  //fin[1] = TFile::Open(("/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalCracks/git"+version+"/version100/model4/gamma/phi_0.440pi/eta20_et60_pu0_IC3_Si2.root").c_str());
  fin[1] = TFile::Open(("/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalCracks/git"+version+"/version100/model4/gamma/eta20_et60_pu0_IC3_Si2.root").c_str());

  if (!fin[0] || !fin[1]) {
    std::cout << " Input files not found. Check paths." << std::endl;
    return 1;
  }

  const unsigned nP = (version.find("V06-")!=version.npos || version.find("V06a")!=version.npos)?62:102;//92
  const unsigned nPuniq = (version.find("V06-")!=version.npos || version.find("V06a")!=version.npos)?62:92;
  const double stepsize = 510./nP;//460./nP;
  double xmin = -170;
  if (linedup) xmin += 0.33+2.5;

  TLatex lat;
  char buf[100];

  TCanvas *mycE[nF];
  mycE[0] = new TCanvas("mycE0","mycE0",1);
  mycE[1] = new TCanvas("mycE1","mycE1",1);

  TCanvas *mycCalib[nF];
  TCanvas *myc[3*nF];
  for (unsigned iF(0); iF<3*nF;++iF){//loop on files
    std::ostringstream lname;
    lname << "mycCalib" << iF;
    if (iF<2) mycCalib[iF] = new TCanvas(lname.str().c_str(),lname.str().c_str(),1500,1000);
    lname.str("");
    lname << "myc" << iF;
    myc[iF] = new TCanvas(lname.str().c_str(),lname.str().c_str(),1500,1000);
  }

  TGraphErrors *gr[nF];
  TGraphErrors *grMean[nF];

  TH1F *hE[nF][nP];
  TH1F *hECalib[nF];
  TH1F *hETot[nF];

  for (unsigned iF(0); iF<nF;++iF){//loop on files
    mycE[iF]->Print(("PlotsCracks/EtotalAll_"+label[iF]+".pdf[").c_str());
    mycCalib[iF]->cd();
    gStyle->SetOptStat("eMRu");
    gStyle->SetOptFit(1111);
    gStyle->SetStatH(0.2);
    //gStyle->SetStatX(0.15);
    fcalib[iF]->cd("Energies");
    
    TTree *treeCalib = (TTree*)gDirectory->Get("Ereso");
    std::ostringstream lcalib;
    lcalib << "Ecalib_" << label[iF];
    hECalib[iF] = new TH1F(lcalib.str().c_str(),";E (Mips); showers",100,20000,28000);
    treeCalib->Draw(("wgtEtotal>>"+lcalib.str()).c_str());
    hECalib[iF]->Fit("gaus");
    TF1 *fit = (TF1*)hECalib[iF]->GetFunction("gaus");
    if (!fit) continue;
    double normalisation = Egen/fit->GetParameter(1);

    double sigmaeff_ref = effSigmaMacro(hECalib[iF]);

    sprintf(buf,"E_{gen} = %3.3f GeV",Egen);
    lat.DrawLatexNDC(0.2,0.8,buf);
    sprintf(buf,"#sigma_{eff} = %3.3f Mips",sigmaeff_ref);
    lat.DrawLatexNDC(0.2,0.7,buf);
    sprintf(buf,"#sigma_{eff} = %3.3f #pm %3.3f GeV",sigmaeff_ref*normalisation,hECalib[iF]->GetRMSError()*normalisation);
    lat.DrawLatexNDC(0.2,0.6,buf);
    sprintf(buf,"#sigma_{eff}/E_{gen} = %3.3f",sigmaeff_ref*normalisation/Egen);
    lat.DrawLatexNDC(0.2,0.5,buf);
    mycCalib[iF]->Update();
    mycCalib[iF]->Print(("PlotsCracks/EtotalCalib_"+label[iF]+".pdf").c_str());

    myc[iF]->Divide(8,8);

    fin[iF]->cd("Energies");
    TTree *tree = (TTree*)gDirectory->Get("Ereso");
    double meandiff[nP];
    double errmeandiff[nP];
    double sigeff[nP];
    double errsig[nP];
    double vtx[nP];
    double errvtx[nP];
    double gaplay[nP];
    double errlay[nP];
    lcalib.str("");
    lcalib << "Etot_" << label[iF];
    hETot[iF] = new TH1F(lcalib.str().c_str(),";E (GeV); showers",6000,0,300);
    for (unsigned iP(0);iP<nP;++iP){//loop on points
      myc[iF]->cd(iP+1);

      std::ostringstream lname,lvar,lcut;
      lname << "Etotal_" << label[iF] << "_step" << iP;

      hE[iF][iP] = new TH1F(lname.str().c_str(),";E (GeV); showers",150,0,300);

      lvar << "wgtEtotal*" << normalisation << ">>" << lname.str();
      lcut << "vtxX>" << xmin+(iP%nPuniq)*stepsize 
	   << " && vtxX<" << xmin+(iP%nPuniq+1)*stepsize;
      errvtx[iP] = stepsize/2.;
      vtx[iP] = xmin+(iP+0.5)*stepsize;
      gaplay[iP] = getGapLayer(vtx[iP]);
      errlay[iP] = 0;
      tree->Draw(lvar.str().c_str(),lcut.str().c_str());

      mycE[iF]->cd();
      hE[iF][iP]->Draw();
      mycE[iF]->Update();
      TPaveStats *stats = (TPaveStats*)(hE[iF][iP]->GetListOfFunctions()->FindObject("stats"));
      stats->SetX1NDC(0.15);
      stats->SetX2NDC(0.51);
      stats->SetY1NDC(0.12);
      stats->SetY2NDC(0.72);

      //double minfit = hE[iF][iP]->GetMean()-2*hE[iF][iP]->GetRMS();
      //double maxfit = hE[iF][iP]->GetMean()+2*hE[iF][iP]->GetRMS();
      hE[iF][iP]->Fit("gaus");//,"RBI","same",minfit,maxfit);
      fit = (TF1*)hE[iF][iP]->GetFunction("gaus");

      lat.DrawLatexNDC(0.2,0.85,lcut.str().c_str());
      if (iF==1){
	lat.DrawLatexNDC(0.2,0.8,"3^{o} tilt");
      }

      sprintf(buf,"#sigma_{eff}=%3.3f #pm %3.3f GeV",effSigmaMacro(hE[iF][iP]),hE[iF][iP]->GetRMSError());
      lat.DrawLatexNDC(0.2,0.75,buf);

      mycE[iF]->Update();
      mycE[iF]->Print(("PlotsCracks/EtotalAll_"+label[iF]+".pdf").c_str());
      double Epeak = hE[iF][iP]->GetMean();
      if (fit) Epeak = fit->GetParameter(1);
      sigeff[iP] = effSigmaMacro(hE[iF][iP])/Epeak;
      errsig[iP] = hE[iF][iP]->GetRMSError()/Epeak;
      meandiff[iP] = (Egen-Epeak)/Egen;
      errmeandiff[iP] = fit? fit->GetParError(1)/Egen : hE[iF][iP]->GetMeanError()/Egen;

      double shift = Egen*1./Epeak;

      lvar.str("");
      lvar << "(wgtEtotal*" << normalisation*shift << ")>>";
      if (iP>0) lvar << "+";
      lvar << lcalib.str();

      tree->Draw(lvar.str().c_str(),lcut.str().c_str());
      //std::cout << "point " << iP << " entries " << hETot[iF]->GetEntries() << " lvar " << lvar.str() << " " << lcut.str()<< " " << hETot[iF]->GetName() << std::endl;
      //hETot[iF]->Draw();
      //return 1;
    }//loop on points

    //return 1;

    myc[iF]->Update();
    myc[iF]->Print(("PlotsCracks/EtotalPerVertexBin_"+label[iF]+".pdf").c_str());


    mycE[iF]->cd();
    hETot[iF]->Draw();
    mycE[iF]->Update();
    TPaveStats* stats = (TPaveStats*)(hETot[iF]->GetListOfFunctions()->FindObject("stats"));
    stats->SetX1NDC(0.15);
    stats->SetX2NDC(0.41);
    stats->SetY1NDC(0.12);
    stats->SetY2NDC(0.42);
    sprintf(buf,"Single #gamma, E_{gen} = %3.1f GeV, #eta=%3.1f",Egen,eta);
    lat.DrawLatexNDC(0.2,0.85,buf);
    if (iF==1){
      lat.DrawLatexNDC(0.2,0.75,"3^{o} tilt");
    }
    double sigefftot = effSigmaMacro(hETot[iF]);
    sprintf(buf,"#sigma_{eff}=%3.3f #pm %3.3f GeV",sigefftot,hETot[iF]->GetRMSError());
    lat.DrawLatexNDC(0.2,0.65,buf);

    int binmin = hETot[iF]->FindBin(Egen-sigmaeff_ref*normalisation);
    int binmax = hETot[iF]->FindBin(Egen+sigmaeff_ref*normalisation);
    double errint = 0;
    double fraction = hETot[iF]->IntegralAndError(binmin,binmax,errint)/hETot[iF]->Integral();
    binmin = hETot[iF]->FindBin(Egen-sigefftot);
    binmax = hETot[iF]->FindBin(Egen+sigefftot);
    double errintref = 0;
    double fractionref = hETot[iF]->IntegralAndError(binmin,binmax,errintref)/hETot[iF]->Integral();

    sprintf(buf,"f_{-#sigma_{ref}}^{#sigma_{ref}}=%3.3f #pm %3.3f (%3.3f #pm %3.3f)",fraction,errint/hETot[iF]->Integral(),fractionref,errintref/hETot[iF]->Integral());
    lat.DrawLatexNDC(0.2,0.5,buf);


    mycE[iF]->Update();
    mycE[iF]->Print(("PlotsCracks/EtotalAllAligned_"+label[iF]+".pdf").c_str());

    gr[iF] = new TGraphErrors(nP,vtx,sigeff,errvtx,errsig);
    //if (linedup) gr[iF] = new TGraphErrors(nP,vtx,sigeff,errvtx,errsig);
    //else gr[iF] = new TGraphErrors(nP,gaplay,sigeff,errlay,errsig);
    gr[iF]->SetMarkerStyle(22);
    gr[iF]->SetMinimum(0);
    gr[iF]->SetMaximum(0.06);
    gr[iF]->SetTitle(";vtx x (mm); #sigma_{eff}/E_{peak}");
    //if (linedup) gr[iF]->SetTitle(";vtx x (mm); #sigma_{eff}/E_{peak}");
    //else gr[iF]->SetTitle(";layer with gap; #sigma_{eff}/E_{peak}");
    myc[2+iF]->cd();
    gStyle->SetOptStat(0);
    gr[iF]->Draw("APE");

    TBox *cr1 = new TBox(-61.7,gr[iF]->GetMinimum(),-51.7,gr[iF]->GetMaximum());
    cr1->SetFillColor(7);
    //cr1->SetFillStyle(3004);
    if (linedup) cr1->Draw();
    TBox *cr2 = new TBox(41.7,gr[iF]->GetMinimum(),51.7,gr[iF]->GetMaximum());
    cr2->SetFillColor(7);
    //cr2->SetFillStyle(3004);
    if (linedup) cr2->Draw();
    TBox *cr3 = new TBox(145,gr[iF]->GetMinimum(),155,gr[iF]->GetMaximum());
    cr3->SetFillColor(7);
    //cr3->SetFillStyle(3004);
    if (linedup) cr3->Draw();

    if (linedup) {
      lat.SetTextColor(9);
      lat.SetTextSize(0.04);
      lat.DrawLatex(-70,0.005,"L10-19 crack");
      lat.DrawLatex(30,0.005,"L20-27 crack");
      lat.DrawLatex(130,0.005,"L0-9 crack");
      lat.SetTextColor(1);
      lat.SetTextSize(0.05);
    }
    else {
      if (version=="V06b-03-05") drawLayerWithGap(lat,random);
      lat.SetTextColor(9);
      lat.SetTextSize(0.03);
      lat.SetTextAngle(90);
      if (version=="V06c-04-06"){
	lat.DrawLatex(230.,0,"Even layers");
	lat.DrawLatex(260.,0,"Odd layers");
      }
      else if (version=="V06d-04-06"){
	lat.DrawLatex(230.,0,"layer%3=0");
	lat.DrawLatex(260.,0,"layer%3=1");
	lat.DrawLatex(290.,0,"layer%3=2");
	lat.DrawLatex(-170.,0,"layer%3=2");
      }
      lat.SetTextAngle(0);
      lat.SetTextColor(1);
      lat.SetTextSize(0.05);
    }

    TLine *ref = new TLine(-170,sigmaeff_ref*normalisation/Egen,290,sigmaeff_ref*normalisation/Egen);
    ref->SetLineColor(2);
    ref->Draw();
    gr[iF]->Draw("PE");

    sprintf(buf,"Single #gamma, E_{gen} = %3.1f GeV, #eta=%3.1f",Egen,eta);
    lat.DrawLatexNDC(0.2,0.8,buf);

    if (iF==1){
      lat.DrawLatexNDC(0.3,0.7,"3^{o} tilt");
    }

    myc[2+iF]->Update();
    myc[2+iF]->Print(("PlotsCracks/SigmaEffvsVtxPos_"+label[iF]+".pdf").c_str());


    grMean[iF] = new TGraphErrors(nP,vtx,meandiff,errvtx,errmeandiff);
    grMean[iF]->SetMarkerStyle(22);
    grMean[iF]->SetMinimum(linedup?-0.04:-0.02);
    grMean[iF]->SetMaximum(linedup?0.16:0.05);
    grMean[iF]->SetTitle(";vtx x (mm); (Egen-Epeak)/E_{gen}");
    myc[4+iF]->cd();
    gStyle->SetOptStat(0);
    grMean[iF]->Draw("APE");

    cr1 = new TBox(-61.7,grMean[iF]->GetMinimum(),-51.7,grMean[iF]->GetMaximum());
    cr1->SetFillColor(7);
    //cr1->SetFillStyle(3004);
    if (linedup) cr1->Draw();
    cr2 = new TBox(41.7,grMean[iF]->GetMinimum(),51.7,grMean[iF]->GetMaximum());
    cr2->SetFillColor(7);
    //cr2->SetFillStyle(3004);
    if (linedup) cr2->Draw();
    cr3 = new TBox(145,grMean[iF]->GetMinimum(),155,grMean[iF]->GetMaximum());
    cr3->SetFillColor(7);
    //cr3->SetFillStyle(3004);
    if (linedup) cr3->Draw();

    if (linedup) {
      lat.SetTextColor(9);
      lat.SetTextSize(0.04);
      lat.DrawLatex(-70,-0.04,"L10-19 crack");
      lat.DrawLatex(30,-0.04,"L20-27 crack");
      lat.DrawLatex(130,-0.04,"L0-9 crack");
      lat.SetTextColor(1);
      lat.SetTextSize(0.05);
    } 
    else { 
      if (version=="V06b-03-05") drawLayerWithGap(lat,random,-0.015);
      lat.SetTextColor(9);
      lat.SetTextSize(0.03);
      lat.SetTextAngle(90);
      if (version=="V06c-04-06"){
	lat.DrawLatex(230.,-0.015,"Even layers");
	lat.DrawLatex(260.,-0.015,"Odd layers");
      }
      else if (version=="V06d-04-06"){
	lat.DrawLatex(230.,-0.015,"layer%3=0");
	lat.DrawLatex(260.,-0.015,"layer%3=1");
	lat.DrawLatex(290.,-0.015,"layer%3=2");
	lat.DrawLatex(-170.,-0.015,"layer%3=2");
      }
      lat.SetTextAngle(0);
      lat.SetTextColor(1);
      lat.SetTextSize(0.05);
    }

    grMean[iF]->Draw("PE");


    sprintf(buf,"Single #gamma, E_{gen} = %3.1f GeV, #eta=%3.1f",Egen,eta);
    lat.DrawLatexNDC(0.2,0.8,buf);

    if (iF==1){
      lat.DrawLatexNDC(0.2,0.7,"3^{o} tilt");
    }


    myc[4+iF]->Update();
    myc[4+iF]->Print(("PlotsCracks/MeanDiffvsVtxPos_"+label[iF]+".pdf").c_str());


    mycE[iF]->Print(("PlotsCracks/EtotalAll_"+label[iF]+".pdf]").c_str());

  }//loop on files



 return 0;
}//main

