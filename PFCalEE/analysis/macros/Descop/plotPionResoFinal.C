#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>

#include "TMath.h"
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

struct FitResult{
  double chi2;
  unsigned ndf;
  double mean;
  double sigma;
  double meanerr;
  double sigmaerr;
};

double E(const unsigned pT, const unsigned eta){
  return pT*cosh(eta/10.);
};

double pT(const unsigned E, const unsigned eta){
  return E/cosh(eta/10.);
};

void drawChi2(TCanvas *myc,TH1F ** p_chi2ndf){
  
  gStyle->SetOptStat("eMRuo");
  gStyle->SetStatH(0.4);
  gStyle->SetStatW(0.4);
  for (unsigned iSR(0); iSR<1;++iSR){
    myc->cd(iSR+1);
    if (p_chi2ndf[iSR]) p_chi2ndf[iSR]->Draw();
  }
  
  myc->Update();
  std::ostringstream lsave;
  lsave << "PionPLOTS/EnergyFitQuality.pdf";
  myc->Print(lsave.str().c_str());
}

double absWeight(const unsigned layer, const bool dedx=true){
  if (dedx==false){
    if (layer == 0) return 0.0378011;//6.285/95.4=0.06588;
    if (layer == 1) return 1;//95.4/95.4=1
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
  }
  else {
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
  }
  return 1;
};


TPad* plot_ratio(TPad *canv, bool up){
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
  if (up){
    pad = new TPad("upper","pad",0, 0.26 ,1 ,1);
    pad->SetBottomMargin(0.05);
    pad->SetTopMargin(0.09);
    pad->Draw();
    pad->cd();
    return pad;
  }
  else {
    pad = new TPad("lower","pad",0, 0   ,1 ,0.26);  
    pad->SetTopMargin(0.05);
    pad->SetBottomMargin(0.24);
    pad->Draw();
    return pad;
  }

};


unsigned fitEnergy(TH1F *hist,
		   TPad *pad,
		   std::string unitStr,
		   FitResult & lres,
		   unsigned isr){
  
  pad->cd();
  //double eMin = hist->GetMean()-5*hist->GetRMS();
  //double eMax = hist->GetMean()+5*hist->GetRMS();
  //hist->GetXaxis()->SetRangeUser(eMin,eMax);
  hist->Draw("PE");


  double nRMSm = 1.5;//isr<1? 1 : 2;
  double nRMSp = 1.5;
  
  TF1 *fitResult = new TF1("fitResult","[0]*TMath::Gaus(x,[1],[2],0)",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
  fitResult->SetParameters(hist->Integral(),
			   hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin()),
			   hist->GetRMS());

  //std::cout << " Initial params: "  << fitResult->GetParameter(0) << " "<< fitResult->GetParameter(1) << " " << fitResult->GetParameter(2)
  //<< std::endl;


  int status = hist->Fit("fitResult","L0QEMI","",
			 fitResult->GetParameter(1)-nRMSm*fitResult->GetParameter(2),
			 fitResult->GetParameter(1)+nRMSp*fitResult->GetParameter(2));
  
  
  //std::cout << " First fit: " << status << " " << fitResult->GetParameter(1) << " " << fitResult->GetParameter(2)
  //<< std::endl;

  // if ((status != 0 && status != 4000) || fitResult->GetChisquare()/fitResult->GetNDF()>5){
  //   //std::cout << " -- Bad fit ! Try again..." << std::endl;
  //   status = hist->Fit("fitResult","L0QEMI","",
  // 		       fitResult->GetParameter(1)-nRMSm*fitResult->GetParameter(2),
  // 		       fitResult->GetParameter(1)+nRMSp*fitResult->GetParameter(2));
    
  //   std::cout << " Second fit: " << status << " " << fitResult->GetParameter(1) << " " << fitResult->GetParameter(2)
  // 	      << std::endl;
  // }
  
  // std::cout << " Final fit: " << fitResult->GetParameter(1) << " " << fitResult->GetParameter(2)
  // 	    << std::endl;
  
  fitResult->SetLineColor(2);
  fitResult->Draw("same");
  
  if (status != 0 && status != 4000) {
    std::cout << " Warning! Fit failed with status " << status << "! Please have a look at the verbose output below...." << std::endl;
    fitResult->SetParameters(hist->Integral(),
			     hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin()),
			     hist->GetRMS());   
    hist->Fit("fitResult","L0QEMI");
    hist->Fit("fitResult","L0EMIV","",
	      fitResult->GetParameter(1)-2*fitResult->GetParameter(2),
	      fitResult->GetParameter(1)+2*fitResult->GetParameter(2));
    //totalE for pu140 is expected to be pathological :/
    //if (isr!=7) return 1;
  }

  char buf[500];
  TLatex lat;
  double latx = hist->GetXaxis()->GetXmin()+(hist->GetXaxis()->GetXmax()-hist->GetXaxis()->GetXmin())/20.;
  double laty = hist->GetMaximum();
  sprintf(buf,"<Efit> = %3.3f +/- %3.3f %s",fitResult->GetParameter(1),fitResult->GetParError(1),unitStr.c_str());
  lat.DrawLatex(latx,laty*0.9,buf);
  sprintf(buf,"RMSfit = %3.3f +/- %3.3f %s",fitResult->GetParameter(2),fitResult->GetParError(2),unitStr.c_str());
  lat.DrawLatex(latx,laty*0.8,buf);
  sprintf(buf,"RMS/meanfit = %3.3f",fitResult->GetParameter(2)/fitResult->GetParameter(1));
  lat.DrawLatex(latx,laty*0.7,buf);
  
  sprintf(buf,"#chi^{2}/N = %3.3f/%d = %3.3f",fitResult->GetChisquare(),fitResult->GetNDF(),fitResult->GetChisquare()/fitResult->GetNDF());
  lat.DrawLatex(latx,laty*0.6,buf);
  
  lres.chi2 = fitResult->GetChisquare();
  lres.ndf = fitResult->GetNDF();
  lres.mean = fitResult->GetParameter(1);
  lres.meanerr = fitResult->GetParError(1);
  lres.sigma = fitResult->GetParameter(2);
  lres.sigmaerr = fitResult->GetParError(2);

  return 0;
};

TPad* plotCalibration(TGraphErrors *gr,TPad *pad,bool doRatio, TGraphErrors *grDelta,std::string unit, double & calib,double & calibErr, double & offset, double & offsetErr,const unsigned eta, const bool dovsE){

  TPad *upper = 0;
  TPad *lower = 0;

  if (!doRatio) pad->cd();
  else {
    pad->Clear();
    upper = plot_ratio(pad, true);
    lower = plot_ratio(pad, false);
    upper->cd();
  }
  gr->GetXaxis()->SetLabelSize(0.0);
  gr->GetXaxis()->SetTitleSize(0.0);
  gr->GetYaxis()->SetLabelSize(0.06);
  gr->GetYaxis()->SetTitleSize(0.06);
  gr->GetXaxis()->SetTitleOffset(7);
  gr->GetYaxis()->SetTitleOffset(0.9);
  
  gr->Draw("ap");

  if (!doRatio){
    if (!dovsE) gr->GetXaxis()->SetTitle("E_{T} (GeV)");
    else gr->GetXaxis()->SetTitle("E (GeV)");
  }
  else gr->GetXaxis()->SetTitle("");
  
  gr->GetYaxis()->SetTitle(("Average energy deposited ("+unit+")").c_str()); 
  char buf[500];
  TF1 *fitFunc=new TF1("calib","[0]+[1]*x",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
  //TF1 *fitFunc=new TF1("calib","[0]+[1]*x+[2]*x*x",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
  fitFunc->SetLineColor(6);

  //if (dovsE) gr->Fit(fitFunc,"RIME","same",0,200);
  //else gr->Fit(fitFunc,"RIME","same",0,pT(200,eta));
  if (dovsE) gr->Fit(fitFunc,"QIME","same");//,0,400);
  else gr->Fit(fitFunc,"QIME","same");
  TLatex lat;
  lat.SetTextColor(6);
  lat.SetTextSize(0.1);
  if (!dovsE) sprintf(buf,"<E> #propto a + b #times E_{T} ");
  else sprintf(buf,"<E> #propto a + b #times E ");
  //else sprintf(buf,"<E> #propto a + b #times E + c #times E^{2} ");
  lat.DrawLatexNDC(0.2,0.85,buf);
  sprintf(buf,"a = %3.3f #pm %3.3f %s",fitFunc->GetParameter(0),fitFunc->GetParError(0),unit.c_str());
  lat.DrawLatexNDC(0.2,0.7,buf);
  sprintf(buf,"b = %3.3f #pm %3.3f %s/GeV",fitFunc->GetParameter(1),fitFunc->GetParError(1),unit.c_str());
  lat.DrawLatexNDC(0.2,0.55,buf);
  //sprintf(buf,"c = %3.3f #pm %3.3f %s/GeV",fitFunc->GetParameter(2),fitFunc->GetParError(2),unit.c_str());
  //lat.DrawLatexNDC(0.2,0.4,buf);
  sprintf(buf,"#chi^{2}/N = %3.3f/%d = %3.3f",fitFunc->GetChisquare(),fitFunc->GetNDF(),fitFunc->GetChisquare()/fitFunc->GetNDF());
  lat.DrawLatexNDC(0.2,0.25,buf);

  calib = fitFunc->GetParameter(1);
  offset = fitFunc->GetParameter(0);
  calibErr = fitFunc->GetParError(1);
  offsetErr = fitFunc->GetParError(0);

  if (doRatio){
    //draw deltaE/E vs E
    lower->cd();
    gPad->SetLogx(0);
    gPad->SetGridx(1);
    gPad->SetGridy(1);

    double loffset = fitFunc->GetParameter(0);
    double lslope = fitFunc->GetParameter(1);
    double range = 0.05;
    if (unit=="GeV") {
      loffset=0;
      lslope=1;
      range = 0.1;
    }

    //fill delta
    for (int ip(0);ip<gr->GetN();++ip){
      double x=0;
      double y=0;
      gr->GetPoint(ip,x,y);
      grDelta->SetPoint(ip,x,((y-loffset)/lslope-x)/x);
      double err = gr->GetErrorY(ip)/lslope*1./x;
      grDelta->SetPointError(ip,0,err);
      //std::cout << "Calib " << ip << " Egen=" << x << " Erec=" << y << " delta=" << ((y-loffset)/lslope-x)/x << std::endl;
    }
    grDelta->SetTitle("");
    grDelta->SetMinimum(-1.*range);
    grDelta->SetMaximum(range);
    grDelta->GetXaxis()->SetLabelSize(0.15);
    grDelta->GetXaxis()->SetTitleSize(0.15);
    grDelta->GetYaxis()->SetLabelSize(0.12);
    grDelta->GetYaxis()->SetTitleSize(0.15);
    grDelta->GetXaxis()->SetTitleOffset(0.5);
    grDelta->GetYaxis()->SetTitleOffset(0.3);
    
    grDelta->Draw("ap");
    //grDelta->GetYaxis()->SetRangeUser(0,grDelta->GetYaxis()->GetXmax());
    if (!dovsE) grDelta->GetXaxis()->SetTitle("E_{T} (GeV)");
    else grDelta->GetXaxis()->SetTitle("E (GeV)");
    grDelta->GetYaxis()->SetTitle("(#Delta E)/E");

    TLine *line = new TLine(grDelta->GetXaxis()->GetXmin(),0,grDelta->GetXaxis()->GetXmax(),0);
    line->SetLineColor(2);//kYellow+4);
    line->Draw();
    

    lower->Update();

  }

  return upper;
};

TPad* plotCalibrationRev(TGraphErrors *gr,TPad *pad,bool doRatio, TGraphErrors *grDelta,std::string unit, 
			 double & calib2,double & calib2Err, 
			 double & calib,double & calibErr, 
			 double & offset, double & offsetErr,const unsigned eta, const bool dovsE){

  TPad *upper = 0;
  TPad *lower = 0;

  if (!doRatio) pad->cd();
  else {
    pad->Clear();
    upper = plot_ratio(pad, true);
    lower = plot_ratio(pad, false);
    upper->cd();
  }
  gr->GetXaxis()->SetLabelSize(0.0);
  gr->GetXaxis()->SetTitleSize(0.0);
  gr->GetYaxis()->SetLabelSize(0.06);
  gr->GetYaxis()->SetTitleSize(0.06);
  gr->GetXaxis()->SetTitleOffset(7);
  gr->GetYaxis()->SetTitleOffset(0.9);
  
  gr->Draw("ap");

  if (!doRatio){
    if (!dovsE) gr->GetYaxis()->SetTitle("E_{T} (GeV)");
    else gr->GetYaxis()->SetTitle("E (GeV)");
  }
  else gr->GetXaxis()->SetTitle("");
  
  gr->GetXaxis()->SetTitle(("Average energy deposited ("+unit+")").c_str()); 
  char buf[500];
  //TF1 *fitFunc=new TF1("calib","[0]+[1]*x",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
  TF1 *fitFunc=new TF1("calib","[0]+[1]*x+[2]*x*x",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
  fitFunc->SetLineColor(6);
  fitFunc->SetParameters(-1.8,0.9,0.0003);

  //if (dovsE) gr->Fit(fitFunc,"RIME","same",0,200);
  //else gr->Fit(fitFunc,"RIME","same",0,pT(200,eta));
  if (dovsE) gr->Fit(fitFunc,"QIME","same");
  else gr->Fit(fitFunc,"QIME","same");
  TLatex lat;
  lat.SetTextColor(6);
  lat.SetTextSize(0.1);
  if (!dovsE) sprintf(buf,"<E> #propto a + b #times E_{T} ");
  //else sprintf(buf,"<E> #propto a + b #times E ");
  else sprintf(buf,"<E> #propto a + b #times E + c #times E^{2} ");
  lat.DrawLatexNDC(0.2,0.85,buf);
  sprintf(buf,"a = %3.3f #pm %3.3f %s",fitFunc->GetParameter(0),fitFunc->GetParError(0),unit.c_str());
  lat.DrawLatexNDC(0.2,0.7,buf);
  sprintf(buf,"b = %3.3f #pm %3.3f %s/GeV",fitFunc->GetParameter(1),fitFunc->GetParError(1),unit.c_str());
  lat.DrawLatexNDC(0.2,0.55,buf);
  sprintf(buf,"c = %3.2fe-5 #pm %3.2fe-5 %s/GeV",fitFunc->GetParameter(2)*100000,fitFunc->GetParError(2)*100000,unit.c_str());
  lat.DrawLatexNDC(0.2,0.4,buf);
  sprintf(buf,"#chi^{2}/N = %3.3f/%d = %3.3f",fitFunc->GetChisquare(),fitFunc->GetNDF(),fitFunc->GetChisquare()/fitFunc->GetNDF());
  lat.DrawLatexNDC(0.2,0.25,buf);

  calib2 = fitFunc->GetParameter(2);
  calib = fitFunc->GetParameter(1);
  offset = fitFunc->GetParameter(0);
  calib2Err = fitFunc->GetParError(2);
  calibErr = fitFunc->GetParError(1);
  offsetErr = fitFunc->GetParError(0);

  if (doRatio){
    //draw deltaE/E vs E
    lower->cd();
    gPad->SetLogx(0);
    gPad->SetGridx(1);
    gPad->SetGridy(1);

    double loffset = fitFunc->GetParameter(0);
    double lslope = fitFunc->GetParameter(1);
    double lslope2 = fitFunc->GetParameter(2);
    double range = 0.05;
    if (unit=="GeV") {
      //loffset=0;
      //lslope=1;
      //lslope2=0;
      range = 0.05;
    }

    //fill delta
    for (int ip(0);ip<gr->GetN();++ip){
      double x=0;
      double y=0;
      gr->GetPoint(ip,x,y);
      grDelta->SetPoint(ip,x,(loffset+lslope*x+lslope2*x*x-y)/y);
      double err = (lslope+2*lslope2*x)*gr->GetErrorX(ip)/y;
      grDelta->SetPointError(ip,0,err);
      //std::cout << "Calib " << ip << " Egen=" << x << " Erec=" << y << " delta=" << ((y-loffset)/lslope-x)/x << std::endl;
    }
    grDelta->SetTitle("");
    grDelta->SetMinimum(-1.*range);
    grDelta->SetMaximum(range);
    grDelta->GetXaxis()->SetLabelSize(0.15);
    grDelta->GetXaxis()->SetTitleSize(0.15);
    grDelta->GetYaxis()->SetLabelSize(0.12);
    grDelta->GetYaxis()->SetTitleSize(0.15);
    grDelta->GetXaxis()->SetTitleOffset(0.5);
    grDelta->GetYaxis()->SetTitleOffset(0.3);
    
    grDelta->Draw("ap");
    //grDelta->GetYaxis()->SetRangeUser(0,grDelta->GetYaxis()->GetXmax());
    if (!dovsE) grDelta->GetXaxis()->SetTitle("E_{T} (GeV)");
    else grDelta->GetXaxis()->SetTitle("E (GeV)");
    grDelta->GetYaxis()->SetTitle("(#Delta E)/E");

    TLine *line = new TLine(grDelta->GetXaxis()->GetXmin(),0,grDelta->GetXaxis()->GetXmax(),0);
    line->SetLineColor(2);//kYellow+4);
    line->Draw();
    

    lower->Update();

  }

  return upper;
};

bool plotResolution(TGraphErrors *gr,TPad *pad,
		    const unsigned ipu,
		    const unsigned eta,
		    const double & stoch0,
		    const double & const0,
		    const double & noise0,
		    double & stoch,double & stochErr, 
		    double & constant, double & constErr,
		    double & noise,double & noiseErr,
		    const bool dovsE){

  pad->cd();
  gr->GetXaxis()->SetLabelSize(0.06);
  gr->GetXaxis()->SetTitleSize(0.06);
  gr->GetYaxis()->SetLabelSize(0.06);
  gr->GetYaxis()->SetTitleSize(0.06);
  gr->GetXaxis()->SetTitleOffset(0.7);
  gr->GetYaxis()->SetTitleOffset(0.8);
  gr->SetMinimum(0);
  gr->SetMaximum(0.4);
  gr->Draw("ap");
  if (!dovsE) gr->GetXaxis()->SetTitle("E_{T} (GeV)");
  else gr->GetXaxis()->SetTitle("E (GeV)");
  gr->GetYaxis()->SetTitle("#sigma/E");

  TF1 *fitFunc =new TF1("reso","sqrt([0]*[0]/x+[1]*[1]+[2]*[2]/(x*x))",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());

  fitFunc->SetParameter(0,stoch0);
  fitFunc->SetParLimits(0,0,1);
  fitFunc->SetParameter(1,const0);
  fitFunc->SetParLimits(1,0,1);
  fitFunc->SetParameter(2,noise0/2.);
  fitFunc->SetParLimits(2,0,noise0);

  if (ipu<3) 
    fitFunc->FixParameter(2,noise0);
  //if (ipu==2) fitFunc->FixParameter(1,const0);
  
  std::cout << " Initial params: "  << fitFunc->GetParameter(0) << " "<< fitFunc->GetParameter(1) << " " << fitFunc->GetParameter(2)
	    << std::endl;

  int status = gr->Fit(fitFunc,"BIME0Q");

  if (fitFunc->GetChisquare()/fitFunc->GetNDF()>20){
    fitFunc->SetParameters(stoch0,const0,noise0);
    if (!dovsE) status = gr->Fit(fitFunc,"BIME0Q","",5,100);
    else status = gr->Fit(fitFunc,"BIME0Q","",10,500);
  }

  fitFunc->SetLineColor(6);
  fitFunc->Draw("same");

  TF1 *fitRef =new TF1("resoRef","sqrt([0]*[0]/x+[1]*[1]+[2]*[2]/(x*x))",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());

  fitRef->SetParameters(0.73,0.029,0.13);
  fitRef->SetLineColor(7);
  fitRef->Draw("same");

  stoch = fitFunc->GetParameter(0);
  stochErr = fitFunc->GetParError(0);
  constant = fitFunc->GetParameter(1);
  constErr = fitFunc->GetParError(1);
  noise = fitFunc->GetParameter(2);
  noiseErr = fitFunc->GetParError(2);

  char buf[500];
  TLatex lat;
  lat.SetTextSize(0.1);
  lat.SetTextColor(6);
  if (!dovsE) sprintf(buf,"#frac{#sigma}{E} #propto #frac{s}{#sqrt{E_{T}}} #oplus c #oplus #frac{n}{E_{T}}");
  else sprintf(buf,"#frac{#sigma}{E} #propto #frac{s}{#sqrt{E}} #oplus c #oplus #frac{n}{E}");

  lat.DrawLatexNDC(0.5,0.85,buf);
  sprintf(buf,"s=%3.3f #pm %3.3f",stoch,stochErr);
  lat.DrawLatexNDC(0.5,0.7,buf);
  sprintf(buf,"c=%3.3f #pm %3.3f",constant,constErr);
  lat.DrawLatexNDC(0.5,0.55,buf);
  sprintf(buf,"n=%3.3f #pm %3.3f",noise,noiseErr);
  lat.DrawLatexNDC(0.5,0.4,buf);
  //sprintf(buf,"status = %d, #chi^{2}/N = %3.1f/%d = %3.1f",status,fitFunc->GetChisquare(),fitFunc->GetNDF(),fitFunc->GetChisquare()/fitFunc->GetNDF());
  sprintf(buf,"#chi^{2}/N = %3.1f/%d = %3.1f",fitFunc->GetChisquare(),fitFunc->GetNDF(),fitFunc->GetChisquare()/fitFunc->GetNDF());
  lat.DrawLatexNDC(0.5,0.25,buf);
  
  if (status != 0 && status != 4000) {
    std::cout << " -- Fit failed with status " << status << std::endl;
    return false;
  }

  return true;
};

double HoverEcor(const double & ee, 
		 const double & fh, 
		 const double & bh,
		 const double & Cglobal,
		 const unsigned model){

  double FHtoEslope = 21.02;
  double FHtoEoffset = 4;
  double BHtoEslope = 8.62;
  double BHtoEoffset = 3.7;
  double ECALslope = 92.83;
  double ECALoffset = 13;
  double BHoverFH = 1.49;
  double HoverE = 5.28;
  double par0 = 30;
  double par1 = 5.31;

  if (model==2812){
    ECALslope = 93.49;
    ECALoffset = -5;
    FHtoEslope = 22.85;
    FHtoEoffset = 0;//-40;
    BHtoEslope = 4.23;
    BHtoEoffset = 0;//-34.5;
    BHoverFH = 1.49;
    par0 = 35;
    par1 = 4.61;
    HoverE = 4.57;
  }
  else if (model==2411){
    ECALslope = 81.64;
    ECALoffset = 2.5;
    FHtoEslope = 22.5;
    FHtoEoffset = 0;//-40;
    BHtoEslope = 5;
    BHtoEoffset = 0;//-34.5;
    BHoverFH = 1.22;
    par0 = 36;
    par1 = 4.15;
    HoverE = 4.07;
  }
  else if (model==1809){
    ECALslope = 66.14;
    ECALoffset = 41;
    FHtoEslope = 18.4;
    FHtoEoffset = 0;//-40;
    BHtoEslope = 5;
    BHtoEoffset = 0;//-34.5;
    BHoverFH = 0.99;
    par0 = 45;
    par1 = 4.15;
    HoverE = 4.1;
  }
   else if (model==28){
    ECALslope = 92.41;
    ECALoffset = 23;
    par0 = 33;
    par1 = 5.38;
  }
  else if (model==24){
    ECALslope = 93.35;
    ECALoffset = 14;
    par0 = 32;//34;
    par1 = 5.18;//5.34;
    FHtoEslope = 21.59;
    FHtoEoffset = -64;
    BHoverFH = 1.53;
    HoverE = 5.15;
  }
  else if (model==18){
    ECALslope = 92.26;
    ECALoffset = 6;
    par0 = 49;//38;
    par1 = 5.85;//5.34;
    FHtoEslope = 19.08;
    FHtoEoffset = -44;
    BHoverFH = 1.30;
    HoverE = 5.80;
  }

  double eecal = (ee-ECALoffset)/ECALslope;
  double fhcal = (fh-FHtoEoffset)/FHtoEslope;
  double bhcal = (bh-BHtoEoffset)/(BHoverFH*BHtoEslope);
  double etot = eecal + (fhcal+bhcal)/HoverE;

  double slope = par0/etot+par1;
  //double par0 = 490;
  //double par1 = 670;
  //double par2 = -5.89;
  //double slope = par2+TMath::Freq((etot-par0)/sqrt(par1*etot));

  return eecal + (fhcal*Cglobal+bhcal)/slope;
  //return eecal + (fhcal*Cglobal+bhcal)/HoverE;

};

int plotPionResoFinal(){//main

  SetTdrStyle();

  bool dovsE = true;

  const unsigned tpmod = 2812;
  std::ostringstream model;
  if (tpmod>0) model << "_tp" << tpmod ;//<< "uni";
  //if (tpmod>0) model << "_ecal" << tpmod;

  const unsigned nRemove = 1;
  std::vector<unsigned> lToRemove;
  unsigned list[nRemove] = {25};//,27,15,1,10,3,18,5,12,7,23,20};
  //const unsigned nGF = 12;
  //unsigned list[nGF] = {1,15,20,22,24,26,28,30,35,40,50};//
  //unsigned list[nGF] = {1,3,4,5,6,7,8,10,12,15,20};
  //unsigned list[nGF] = {1,2,3,5,10,15,20,25,30,35,40,50};

  const unsigned nPu = 3;//4;
  unsigned pu[nPu] = {0,0,0};//,140,200};

  const unsigned nS = 1;
  std::string scenario[nS] = {
    "pi-/"
  };

  unsigned eta=25;
  //const unsigned neta = 7;
  //unsigned eta[neta]={17,19,21,23,25,27,29};

  double etaval=2.5;

  const unsigned nEvtMin = 150;


  const unsigned nV = 1;
  TString version[nV] = {"33"};//,"0"};
  
  const unsigned nLayers = 54;

  const unsigned nSR = 1;
  double fitQual[nSR];
  double srval[nSR];
  double srerr[nSR];
  fitQual[0] = 50;
  for (unsigned iSR(0); iSR<nSR;++iSR){
    srval[iSR] = iSR*1.;
    srerr[iSR] = 0.;
    if (iSR>0) fitQual[iSR] = 30;
  }


  double calib2[nPu][nLayers][nSR];
  double calib2Err[nPu][nLayers][nSR];
  double calib[nPu][nLayers][nSR];
  double calibErr[nPu][nLayers][nSR];
  double offset[nPu][nLayers][nSR];
  double offsetErr[nPu][nLayers][nSR];

  double sigmaStoch[nPu][nSR][nLayers];
  double sigmaStochErr[nPu][nSR][nLayers];
  double sigmaConst[nPu][nSR][nLayers];
  double sigmaConstErr[nPu][nSR][nLayers];
  double sigmaNoise[nPu][nSR][nLayers];
  double sigmaNoiseErr[nPu][nSR][nLayers];
  
  std::ostringstream saveName;

  unsigned genEnAll[]={3,5,10,30,50,70,100,200};
  //unsigned genEnAll[]={7,10,20,30,40};
  const unsigned nGenEnAll=sizeof(genEnAll)/sizeof(unsigned);


  //canvas so they are created only once
  TCanvas *mycE[nGenEnAll];
  for (unsigned iE(0); iE<nGenEnAll;++iE){
    std::ostringstream lName;
    lName << "mycE" << genEnAll[iE];
    mycE[iE] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1500,1000);
    //mycE[iE]->Divide(2,1);
  }

  TCanvas *mycEtot = new TCanvas("mycEtot","mycEtot",1500,1000);
  mycEtot->Divide(4,2);
  TCanvas *mycCalibLayer = new TCanvas("mycCalibLayer","mycCalibLayer",1500,1000);
  TCanvas *mycCalib2Layer = new TCanvas("mycCalib2Layer","mycCalib2Layer",1500,1000);
  TCanvas *mycOffsetLayer = new TCanvas("mycOffsetLayer","mycOffsetLayer",1500,1000);
  
  const unsigned nCanvas = 4;  
  TCanvas *myc[nCanvas];
  for (unsigned iC(0);iC<nCanvas;++iC){
    std::ostringstream lName;
    lName << "myc" << iC;
    myc[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1500,1000);
    //myc[iC]->Divide(1,2);
  }

  TCanvas *mycReso = new TCanvas("mycReso","mycReso",750,500);
  TCanvas *mycR = new TCanvas("mycR","Sampling",1500,1000);
  TCanvas *mycC = new TCanvas("mycC","Constant",1500,1000);
  TCanvas *mycN = new TCanvas("mycN","Noise",1500,1000);

  TH1F *p_chi2ndf[nSR];
  for (unsigned iSR(0); iSR<nSR;++iSR){//loop on signal region
    std::ostringstream label;
    label << "p_chi2ndf_SR" << iSR;
    p_chi2ndf[iSR] = new TH1F(label.str().c_str(),";#chi^{2}/N;entries",500,0,50);
    p_chi2ndf[iSR]->StatOverflows();
  }


  for (unsigned ic(0);ic<nRemove;++ic){//loop on intercalib
    unsigned minlayer = 0;
    double minval = 1000;
    TString pSuffix = "_Final"+model.str()+"_Remove";
    pSuffix += ic;

    TFile *fcalib;
    std::ostringstream label;
    label << "PionPLOTS/CalibReso";
    if (dovsE) label << "_vsE";
    label << "_Remove" << ic;
    label << "_nofit.root";
    fcalib = TFile::Open(label.str().c_str(),"RECREATE");
    
    for (unsigned iV(0); iV<nV;++iV){//loop on versions
      for (unsigned iS(0); iS<nS;++iS){//loop on scenarios
	
	TString plotDir = "/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalDescop/gitV05-02-04/version"+version[iV]+"/"+scenario[iS]+"/";
	TTree *ltree[nPu][nGenEnAll];
	TGraphErrors *resoRecoFit[nPu][nLayers][nSR];
	
	for (unsigned ipu(0); ipu<nPu; ++ipu){//loop on pu
	  unsigned puOption = pu[ipu];
	  
	  std::string unit = "MIPS";
	  if (ipu>0) unit = "GeV";
	  
	  //identify valid energy values
	  bool skip[nGenEnAll];
	  unsigned nValid = 0;
	  TFile *inputFile[nGenEnAll];
	  for (unsigned iE(0); iE<nGenEnAll; ++iE){
	    ltree[ipu][iE] = 0;
	    skip[iE] = false;
	    inputFile[iE] = 0;
	    std::ostringstream linputStr;
	    linputStr << plotDir ;
	    //linputStr << "eta" << eta << "_et" << genEnAll[iE] << "_pu" << pu[ipu] << "_IC3";
	    //linputStr << "_nofit.root";
	    linputStr << "pion_gc3-50_reso" << model.str() << "_et" << genEnAll[iE] << ".root";
	    inputFile[iE] = TFile::Open(linputStr.str().c_str());
	    if (!inputFile[iE]) {
	      std::cout << " -- Error, input file " << linputStr.str() << " cannot be opened. Skipping..." << std::endl;
	      //	    return 1;
	      skip[iE] = true;
	    }
	    else {
	      //inputFile[iE]->cd("Energies");
	      ltree[ipu][iE] = (TTree*)gDirectory->Get("Ereso");
	      
	      if (!ltree[ipu][iE]){
		std::cout << " -- File " << inputFile[iE]->GetName() << " sucessfully opened but tree Ereso not found! Skipping." << std::endl;
		skip[iE] = true;
	      } else { 
	      std::cout << " -- File " << inputFile[iE]->GetName() << " sucessfully opened and tree found." << std::endl;
	      nValid++;
	      }
	    }
	  }
	  
	  unsigned newidx = 0;
	  const unsigned nGenEn = nValid;
	  unsigned genEn[nGenEn];
	  unsigned oldIdx[nGenEn];
	  for (unsigned iE(0); iE<nGenEnAll; ++iE){	  
	    if (!skip[iE]) {
	      genEn[newidx]=genEnAll[iE];
	      oldIdx[newidx]=iE;
	      newidx++;
	    }
	  }

	  minval = 1000;
	  unsigned iL=0;
	  //for (unsigned iL(0);iL<nLayers;++iL){//loop on layers
	    
	    TH1F *p_Ereco[nGenEn][nSR];
	    TGraphErrors *calibRecoFit[nSR];
	    TGraphErrors *calibRecoDelta[nSR];
	    
	    //draw calibration curves
	    for (unsigned iSR(0); iSR<nSR;++iSR){//loop on signal region
	      TString srStr = "";
	      srStr += iSR;
	      srStr += "lay";
	      srStr += iL;
	      srStr += "pu";
	      srStr += ipu;
	      calibRecoFit[iSR] = new TGraphErrors();
	      calibRecoFit[iSR]->SetName("calibRecoFit"+srStr);
	      calibRecoFit[iSR]->SetTitle("");
	      calibRecoFit[iSR]->SetMarkerStyle(20);
	      calibRecoFit[iSR]->SetMarkerColor(1);
	      calibRecoFit[iSR]->SetLineColor(1);
	      calibRecoDelta[iSR] = (TGraphErrors *) calibRecoFit[iSR]->Clone("calibRecoDelta"+srStr);
	      resoRecoFit[ipu][iL][iSR] = (TGraphErrors *) calibRecoFit[iSR]->Clone("resoRecoFit"+srStr);
	    }
	    
	    for (unsigned iE(0); iE<nGenEn; ++iE){
	      
	      std::cout << "- Processing energy : " << genEn[iE] 
			<< std::endl;
	      
	      //inputFile[oldIdx[iE]]->cd("Energies");
	      
	      std::cout << " -- Tree entries for eta=" << eta << " pu=" << pu[ipu] << " : " << ltree[ipu][oldIdx[iE]]->GetEntries() << std::endl;
	      
	      for (unsigned iSR(0); iSR<nSR;++iSR){//loop on signal region
		//std::cout << " --Processing signal region: " << iSR << std::endl;

		mycE[iE]->cd();//iSR+1);
	      gStyle->SetOptStat(0);
	      gStyle->SetOptFit(0);

	      std::ostringstream lName;
	      lName.str("");
	      if (ipu>1) lName << offset[1][iL][iSR] << "+" << calib[1][iL][iSR] << "*";
	      //if (ipu>0) lName << "(";
	      if (genEn[iE]<60) 
	      lName << "HoverEcor(EECAL,EFHCAL,EBHCAL,globalC_15," << tpmod << ")";
	      else if (genEn[iE]<300) 
	      lName << "HoverEcor(EECAL,EFHCAL,EBHCAL,globalC_15," << tpmod << ")";
	      else lName << "HoverEcor(EECAL,EFHCAL,EBHCAL,globalC_15)";
	      //if (ipu>0) lName << " - " << offset[0][iL][iSR] << ")/" << calib[0][iL][iSR];

	      if (ipu>1) {
		lName << "+" << calib2[1][iL][iSR] << "*pow(";
		if (genEn[iE]<60) 
		  lName << "HoverEcor(EECAL,EFHCAL,EBHCAL,globalC_15," << tpmod << ")";
		else if (genEn[iE]<300) 
		  lName << "HoverEcor(EECAL,EFHCAL,EBHCAL,globalC_15," << tpmod << ")";
		else lName << "HoverEcor(EECAL,EFHCAL,EBHCAL,globalC_15)";
		//lName << " - " << offset[0][iL][iSR] << ")/" << calib[0][iL][iSR];
		lName << ",2)";
	      }

	      //std::cout << lName.str() << std::endl;
	      ltree[ipu][oldIdx[iE]]->Draw(lName.str().c_str(),"","");
	      lName.str("");
	      lName << "energy" << genEn[iE] << "_SR" << iSR ;
	      p_Ereco[iE][iSR] = (TH1F*)(gPad->GetPrimitive("htemp"))->Clone(lName.str().c_str()); // 1D
	      //if (iSR==7) p_Ereco[iE][iSR] = (TH1F*)gDirectory->Get("p_wgtEtotal");
	      if (!p_Ereco[iE][iSR]){
		std::cout << " -- ERROR, pointer for histogram " << lName.str() << " is null." << std::endl;
		return 1;
	      }
	      //std::cout << " --- Reco E = entries " << p_Ereco[iE][iSR]->GetEntries() 
	      //<< " mean " << p_Ereco[iE][iSR]->GetMean() 
	      //<< " rms " << p_Ereco[iE][iSR]->GetRMS() 
	      //<< " overflows " << p_Ereco[iE][iSR]->GetBinContent(p_Ereco[iE][iSR]->GetNbinsX()+1)
	      //<< std::endl;
	      
	      
	      p_Ereco[iE][iSR]->SetTitle((";E ("+unit+");events").c_str());

	      //take min 20 bins
	      //if(p_Ereco[iE][iSR]->GetNbinsX()>40) p_Ereco[iE][iSR]->Rebin(2);

	      //skip data with too little stat
	      if (p_Ereco[iE][iSR]->GetEntries()<nEvtMin) {
		gPad->Clear();
		continue;
	      }

	      TPad *lpad = (TPad*)(mycE[iE]->cd());//iSR+1));
	      FitResult lres;
	      if (fitEnergy(p_Ereco[iE][iSR],lpad,unit,lres,iSR)!=0) return 1;
	      lpad->cd();
	      char buf[500];
	      sprintf(buf,"#pi^{-} E_{T}=%d GeV + PU %d",genEn[iE],pu[ipu]);
	      TLatex lat;
	      lat.SetTextSize(0.05);
	      lat.DrawLatexNDC(0.25,0.965,buf);
	      sprintf(buf,"#eta=%3.1f",etaval);
	      lat.SetTextSize(0.06);
	      lat.DrawLatexNDC(0.15,0.87,buf);
	      //if (iSR==3) 
		lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");
	      p_chi2ndf[iSR]->Fill(lres.chi2/lres.ndf);

	      //filter out bad points
	      if (lres.chi2/lres.ndf > fitQual[iSR]) {
		std::cout << " --- INFO! Point Egen=" 
			  << genEn[iE] 
			  << " eta=" << etaval
			  << " pu=" << pu[ipu]
			  << " skipped, chi2/ndf = "
			  << lres.chi2/lres.ndf
			  << std::endl;
		continue;
	      }

	      Int_t np=calibRecoFit[iSR]->GetN();
	      if (ipu<1){
		if (!dovsE) calibRecoFit[iSR]->SetPoint(np,genEn[iE],lres.mean);
		else calibRecoFit[iSR]->SetPoint(np,E(genEn[iE],eta),lres.mean);
		calibRecoFit[iSR]->SetPointError(np,0.0,lres.meanerr);
	      } else {
		if (!dovsE) calibRecoFit[iSR]->SetPoint(np,lres.mean,genEn[iE]);
		else calibRecoFit[iSR]->SetPoint(np,lres.mean,E(genEn[iE],eta));
		calibRecoFit[iSR]->SetPointError(np,lres.meanerr,0.0);
	      }
	      double reso = fabs(lres.sigma/lres.mean);
	      if (!dovsE) resoRecoFit[ipu][iL][iSR]->SetPoint(np,genEn[iE],reso);
	      else resoRecoFit[ipu][iL][iSR]->SetPoint(np,E(genEn[iE],eta),reso);
	      double errFit = reso*sqrt(pow(lres.sigmaerr/lres.sigma,2)+pow(lres.meanerr/lres.mean,2));
	      resoRecoFit[ipu][iL][iSR]->SetPointError(np,0,errFit);

	    }//loop on SR

	    saveName.str("");
	    saveName << plotDir << "/Ereco_lay" << iL << "_pu" << puOption;
	    if (ipu==0) saveName << "raw";
	    saveName << "_E" << genEn[iE] ;
	    mycE[iE]->Update();
	    mycE[iE]->Print((saveName.str().c_str()+pSuffix)+".pdf");
	    
	  }//loop on energies
	 
	    //drawChi2(myc[2],p_chi2ndf);


	  //plot and fit calib
	  if (ipu==2){
	    for (unsigned iSR(0); iSR<nSR;++iSR){//loop on signal region
	      if (pu[ipu]!=0 && iSR==nSR-1) continue;
	      TPad *lpad = (TPad*)(myc[ipu]->cd());//iSR+1));
	      TPad *upper = plotCalibration(calibRecoFit[iSR],lpad,
					    true,calibRecoDelta[iSR],
					    unit,
					    calib[ipu][iL][iSR],
					    calibErr[ipu][iL][iSR],
					    offset[ipu][iL][iSR],
					    offsetErr[ipu][iL][iSR],
					    eta,dovsE);
	      upper->cd();
	      char buf[500];
	      sprintf(buf,"#pi^{-} #eta=%3.1f + PU %d",etaval,pu[ipu]);
	      TLatex lat;
	      lat.SetTextSize(0.1);
	      lat.DrawLatexNDC(0.7,0.15,buf);
	      //sprintf(buf,"SR %d",iSR);
	      //lat.DrawLatexNDC(0.7,0.3,buf);
	      lat.SetTextSize(0.06);
	      if (iSR==0) lat.DrawLatexNDC(0.01,0.94,"HGCAL G4 standalone");
	      std::ostringstream lsave;
	      lsave.str("");
	      lsave << "SR" << iSR;
	      lsave << "_remove" << iL;
	      fcalib->mkdir(lsave.str().c_str());
	      fcalib->cd(lsave.str().c_str());
	      calibRecoFit[iSR]->Write();
	      calibRecoDelta[iSR]->Write();
	      //resoRecoFit[ipu][iL][iSR]->Write();

	    }

	    myc[ipu]->Update();
	    std::ostringstream lsave;
	    lsave.str("");
	    lsave << plotDir << "/";
	    if (ipu==0) lsave << "CalibMipToGeV";
	    else lsave << "Calib";
	    lsave << "_lay" << iL << "_pu" << puOption << pSuffix;
	    if (dovsE) lsave << "_vsE";
	    myc[ipu]->Print((lsave.str()+".pdf").c_str());
	    }
	  else if (ipu==1){
	    for (unsigned iSR(0); iSR<nSR;++iSR){//loop on signal region
	      if (pu[ipu]!=0 && iSR==nSR-1) continue;
	      TPad *lpad = (TPad*)(myc[3]->cd());//iSR+1));
	      TPad *upper = plotCalibrationRev(calibRecoFit[iSR],lpad,
					       true,calibRecoDelta[iSR],
					       unit,
					       calib2[ipu][iL][iSR],
					       calib2Err[ipu][iL][iSR],
					       calib[ipu][iL][iSR],
					       calibErr[ipu][iL][iSR],
					       offset[ipu][iL][iSR],
					       offsetErr[ipu][iL][iSR],
					       eta,dovsE);
	      upper->cd();
	      char buf[500];
	      sprintf(buf,"#pi^{-} #eta=%3.1f + PU %d",etaval,pu[ipu]);
	      TLatex lat;
	      lat.SetTextSize(0.1);
	      lat.DrawLatexNDC(0.7,0.15,buf);
	      //sprintf(buf,"SR %d",iSR);
	      //lat.DrawLatexNDC(0.7,0.3,buf);
	      lat.SetTextSize(0.06);
	      if (iSR==0) lat.DrawLatexNDC(0.01,0.94,"HGCAL G4 standalone");
	      std::ostringstream lsave;
	      lsave.str("");
	      lsave << "SR" << iSR;
	      lsave << "_remove" << iL;
	      fcalib->mkdir(lsave.str().c_str());
	      fcalib->cd(lsave.str().c_str());
	      calibRecoFit[iSR]->Write();
	      calibRecoDelta[iSR]->Write();
	      //resoRecoFit[ipu][iL][iSR]->Write();

	    }

	    myc[3]->Update();
	    std::ostringstream lsave;
	    lsave.str("");
	    lsave << plotDir << "/";
	    lsave << "CalibCor";
	    lsave << "_lay" << iL << "_pu" << puOption << pSuffix;
	    if (dovsE) lsave << "_vsE";
	    myc[3]->Print((lsave.str()+".pdf").c_str());
	    }
	  if (ipu==2) {
	    //plot reso
	    for (unsigned iSR(0); iSR<nSR;++iSR){//loop on signal region
	      if (pu[ipu]!=0 && iSR==nSR-1) continue;
	      TPad *lpad = (TPad*)(myc[1]->cd());//iSR+1));
	      
	      double stoch0 = pu[ipu]==0? (dovsE?0.25 : 0.14) : sigmaStoch[1][iSR][iL];
	      double const0 = pu[ipu]==0? 0.01 : sigmaConst[1][iSR][iL];
	      
	      //limit range to get more realistic RMS ?
	      //if (pu[ipu]!=0 && iSR<(nSR-1)) p_sigma[iSR]->GetXaxis()->SetRangeUser(0,5);
	      double noise0 = 0.13;//pu[ipu]==0? 0 : p_sigma[iSR]->GetRMS();
	      
	      
	      bool success = plotResolution(resoRecoFit[ipu][iL][iSR],lpad,
					    ipu,eta,
					    stoch0,const0,noise0,
					    sigmaStoch[ipu][iSR][iL],
					    sigmaStochErr[ipu][iSR][iL],
					    sigmaConst[ipu][iSR][iL],
					    sigmaConstErr[ipu][iSR][iL],
					    sigmaNoise[ipu][iSR][iL],
					    sigmaNoiseErr[ipu][iSR][iL],
					    dovsE);
	      lpad->cd();
	      char buf[500];
	      sprintf(buf,"#pi^{-} #eta=%3.1f + PU %d",etaval,pu[ipu]);
	      TLatex lat;
	      lat.SetTextSize(0.06);
	      lat.DrawLatexNDC(0.15,0.85,buf);
	      //sprintf(buf,"SR %d",iSR);
	      //lat.DrawLatexNDC(0.15,0.7,buf);
	      //if (iSR==4) 
	      lat.SetTextSize(0.06);
		lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");
	      
	      if (!success) {
		return 1;	    
	      }
	      
	      std::ostringstream lsave;
	      lsave.str("");
	      lsave << "SR" << iSR;
	      lsave << "_remove" << iL;
	      //fcalib->mkdir(lsave.str().c_str());
	      fcalib->cd(lsave.str().c_str());
	      //calibRecoFit[iSR]->Write();
	      //calibRecoDelta[iSR]->Write();
	      resoRecoFit[ipu][iL][iSR]->Write();
	      
	    }//loop on SR
	    if (ipu==1 && iL>0 && iL<nLayers-2 && sigmaConst[ipu][4][iL]<minval){
	      bool ignore = false;
	      for (unsigned r(0); r<lToRemove.size(); r++){
		if (abs(static_cast<int>(iL)-static_cast<int>(lToRemove[r]))<2) ignore = true;
	      }
	      if (!ignore){
		minlayer = iL;
		minval = sigmaConst[ipu][4][iL];
	      }
	    }
	    myc[1]->Update();
	    std::ostringstream lsave;
	    lsave.str("");
	    lsave << plotDir << "/";
	    if (ipu==0) lsave << "ResoRaw";
	    else lsave << "Reso";
	    lsave << "_lay" << iL << "_pu" << puOption << pSuffix;
	    if (dovsE) lsave << "_vsE";
	    myc[1]->Print((lsave.str()+".pdf").c_str());

	    }
	  
	  //}//loop on layers
	  
	}//loop on pu

      if (nGenEnAll==1) continue;
      //return 1;

      TLatex lat;
      char buf[500];

      
      //plot and fit calib parameters vs layer removed
      TGraphErrors *grSlope[nPu][nSR];
      TGraphErrors *grSlope2[nPu][nSR];
      TGraphErrors *grOffset[nPu][nSR];

      gStyle->SetOptStat(0);
      gStyle->SetOptFit(0);//1111);
      //gStyle->SetStatX(0.64);//top-right corner
      //gStyle->SetStatY(0.85);
      //gStyle->SetStatW(0.25);
      //gStyle->SetStatH(0.25);
      //gStyle->SetStatColor(0);

      for (unsigned ipu(0); ipu<1; ++ipu){//loop on pu
	mycCalib2Layer->Clear();
	mycCalibLayer->Clear();
	mycOffsetLayer->Clear();
	//mycCalib2Layer->Divide(1,2);
	//mycCalibLayer->Divide(1,2);
	//mycOffsetLayer->Divide(1,2);

	for (unsigned iSR(0); iSR<nSR;++iSR){//loop on signal region
	  if (pu[ipu]!=0 && iSR==nSR-1) continue;
	  grSlope[ipu][iSR] = new TGraphErrors();
	  grSlope2[ipu][iSR] = new TGraphErrors();
	  std::ostringstream llabel;
	  llabel << "grSlope_" << ipu << "_SR" << iSR;
	  grSlope[ipu][iSR]->SetName(llabel.str().c_str());
	  llabel.str("");
	  llabel << ";layer removed;Calib slope (" ;
	  if (ipu==0) llabel << "MIPs/GeV)";
	  else llabel << "GeV/GeV)";
	  grSlope[ipu][iSR]->SetTitle(llabel.str().c_str());
	  grSlope[ipu][iSR]->SetMarkerStyle(20+iSR);
	  grSlope[ipu][iSR]->SetMarkerColor(iSR>3?2+iSR:1+iSR);
	  grSlope[ipu][iSR]->SetLineColor(iSR>3?2+iSR:1+iSR);
	  llabel.str("");
	  llabel << "grSlope2_" << ipu << "_SR" << iSR;
	  grSlope2[ipu][iSR]->SetName(llabel.str().c_str());
	  llabel.str("");
	  llabel << ";layer removed;Calib slope^{2} (" ;
	  if (ipu==0) llabel << "MIPs/GeV^2)";
	  else llabel << "GeV/GeV^2)";
	  grSlope2[ipu][iSR]->SetTitle(llabel.str().c_str());
	  grSlope2[ipu][iSR]->SetMarkerStyle(20+iSR);
	  grSlope2[ipu][iSR]->SetMarkerColor(iSR>3?2+iSR:1+iSR);
	  grSlope2[ipu][iSR]->SetLineColor(iSR>3?2+iSR:1+iSR);
	  llabel.str("");
	  llabel << "grOffset_" << ipu << "_SR" << iSR;
	  grOffset[ipu][iSR] = (TGraphErrors *)grSlope[ipu][iSR]->Clone(llabel.str().c_str());
	  llabel.str("");
	  llabel << ";layer removed;Calib offset (" ;
	  if (ipu==0) llabel << "MIPs)";
	  else llabel << "GeV)";
	  grOffset[ipu][iSR]->SetTitle(llabel.str().c_str());

	  for (unsigned iL(0);iL<1;++iL){//loop on layers
	    grSlope2[ipu][iSR]->SetPoint(iL,iL,calib2[ipu][iL][iSR]);
	    grSlope2[ipu][iSR]->SetPointError(iL,0.0,calib2Err[ipu][iL][iSR]);

	    grSlope[ipu][iSR]->SetPoint(iL,iL,calib[ipu][iL][iSR]);
	    grSlope[ipu][iSR]->SetPointError(iL,0.0,calibErr[ipu][iL][iSR]);

	    grOffset[ipu][iSR]->SetPoint(iL,iL,offset[ipu][iL][iSR]);
	    grOffset[ipu][iSR]->SetPointError(iL,0.0,offsetErr[ipu][iL][iSR]);
	  }//loop on layers
	  mycCalibLayer->cd(iSR+1);
	  grSlope[ipu][iSR]->Draw();
	  //grSlope[ipu][iSR]->Fit("pol2","","same");
	     
	  lat.SetTextSize(0.06);
	  sprintf(buf,"pu %d",pu[ipu]);
	  lat.DrawLatexNDC(0.2,0.87,buf);
	  //if (iSR==4) 
	    lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");

	  mycCalib2Layer->cd(iSR+1);
	  grSlope2[ipu][iSR]->Draw();
	  //grSlope[ipu][iSR]->Fit("pol2","","same");
	     
	  lat.SetTextSize(0.06);
	  sprintf(buf,"pu %d",pu[ipu]);
	  lat.DrawLatexNDC(0.2,0.87,buf);
	  //if (iSR==4) 
	    lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");

	  mycOffsetLayer->cd(iSR+1);
	  grOffset[ipu][iSR]->Draw();
	  //grOffset[ipu][iSR]->Fit("pol2","","same");

	  lat.DrawLatexNDC(0.2,0.87,buf);
	  //if (iSR==4) 
	    lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");

	  std::ostringstream lsave;
	  lsave.str("");
	  lsave << "SR" << iSR;
	  fcalib->mkdir(lsave.str().c_str());
	  fcalib->cd(lsave.str().c_str());
	  grOffset[ipu][iSR]->Write();
	  grSlope[ipu][iSR]->Write();
	  grSlope2[ipu][iSR]->Write();
	}//loop on SR
	mycCalibLayer->Update();
	std::ostringstream lsave;
	lsave << plotDir << "/CalibSlopevsLayer_pu" << ipu;
	if (ipu==0) lsave << "raw";
	if (dovsE) lsave << "_vsE";
	lsave  << pSuffix << ".pdf";
	mycCalibLayer->Print(lsave.str().c_str());
	mycCalib2Layer->Update();
	lsave.str("");
	lsave << plotDir << "/CalibSlope2vsLayer_pu" << ipu;
	if (ipu==0) lsave << "raw";
	if (dovsE) lsave << "_vsE";
	lsave  << pSuffix << ".pdf";
	mycCalib2Layer->Print(lsave.str().c_str());
	lsave.str("");
	lsave << plotDir << "/CalibOffsetvsLayer_pu" << ipu;
	if (ipu==0) lsave << "raw";
	if (dovsE) lsave << "_vsE";
	lsave << pSuffix << ".pdf";
	mycOffsetLayer->Update();
	mycOffsetLayer->Print(lsave.str().c_str());

      }//loop on pu

      TGraph *grDummy = new TGraph();
      grDummy->SetName("grDummy");
      if (dovsE){
	grDummy->SetPoint(0,E(genEnAll[0],eta),0.1);
	grDummy->SetPoint(1,E(genEnAll[nGenEnAll-1],eta),0.1);
      } else {
	grDummy->SetPoint(0,genEnAll[0],0.1);
	grDummy->SetPoint(1,genEnAll[nGenEnAll-1],0.1);
      }
      grDummy->SetLineColor(10);
      grDummy->SetMarkerColor(10);
      grDummy->SetMinimum(0);
      if (!dovsE) grDummy->GetXaxis()->SetTitle("E_{T} (GeV)");
      else grDummy->GetXaxis()->SetTitle("E (GeV)");
      grDummy->GetYaxis()->SetTitle("#sigma/E");
      grDummy->GetYaxis()->SetTitleOffset(1.2);
      gStyle->SetOptStat(0);
      gStyle->SetOptFit(0);
      for (unsigned ipu(1); ipu<nPu; ++ipu){//loop on pu
	if (pu[ipu]==0) grDummy->SetMaximum(0.4);
	else grDummy->SetMaximum(0.16);
	for (unsigned iSR(0); iSR<nSR;++iSR){//loop on signal region
	  if (pu[ipu]!=0 && iSR==nSR-1) continue;
	  
	  mycReso->cd();
	  gPad->SetLogx(1);
	  gPad->SetGridx(1);
	  gPad->SetGridy(1);
	  TLegend *leglay =  new TLegend(0.8,0.5,0.94,0.94);
	  leglay->SetFillColor(10);
	  grDummy->Draw("AP");
	  for (unsigned iL(0);iL<1;++iL){//loop on lay
	    resoRecoFit[ipu][iL][iSR]->SetLineColor(iL%9+1);
	    resoRecoFit[ipu][iL][iSR]->SetMarkerColor(iL%9+1);
	    resoRecoFit[ipu][iL][iSR]->SetMarkerStyle(iL/9+20);
	    resoRecoFit[ipu][iL][iSR]->Draw("PLsame");
	    label.str("");
	    //label << "lay = " << iL;
	    label << "C_{global} = " << list[iL];
	    leglay->AddEntry(resoRecoFit[ipu][iL][iSR],label.str().c_str(),"P");
	  }//loop on lay
	  //add ref pu=0 curve
	  if (pu[ipu]!=0) {
	    resoRecoFit[1][0][iSR]->SetLineWidth(2);
	    resoRecoFit[1][0][iSR]->SetLineStyle(2);
	    resoRecoFit[1][0][iSR]->Draw("Lsame");
	    leglay->AddEntry(resoRecoFit[1][0][iSR],"Ref pu=0","L");
	  }
	  leglay->Draw("same");
	  sprintf(buf,"#pi^{-} + PU %d",pu[ipu]);
	  //lat.SetTextSize(0.1);
	  lat.DrawLatexNDC(0.45,0.85,buf);
	  //sprintf(buf,"SR %d",iSR);
	  //lat.DrawLatexNDC(0.75,0.5,buf);
	  lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");
	  mycReso->Update();
	  std::ostringstream lsave;
	  lsave << plotDir << "/Resolution_pu" << pu[ipu] << "_SR" << iSR;
	  if (dovsE) lsave << "_vsE";
	  lsave << pSuffix;
	  mycReso->Print((lsave.str()+".pdf").c_str());
	  mycReso->Print((lsave.str()+".C").c_str());
	  
	}
      }


      TLegend *leg = new TLegend(0.85,0.7,0.94,0.94);
      leg->SetFillColor(10);

      /*mycR->Divide(2,1);
      mycC->Divide(2,1);
      mycN->Divide(2,1);
      TPad *mypad[3][neta];
      TPad *left[3];
      TPad *right[3];
      left[0] = (TPad*)mycR->cd(1);
      right[0] = (TPad*)mycR->cd(2);
      left[1] = (TPad*)mycC->cd(1);
      right[1] = (TPad*)mycC->cd(2);
      left[2] = (TPad*)mycN->cd(1);
      right[2] = (TPad*)mycN->cd(2);
      for (unsigned iC(0);iC<3;++iC){
	left[iC]->Divide(1,4);
	right[iC]->Divide(1,3);
	for (unsigned ieta=0; ieta<neta;++ieta){//loop on pt values
	  if (ieta<4) mypad[iC][ieta] = (TPad*)left[iC]->GetPad(ieta+1);
	  else mypad[iC][ieta] = (TPad*)right[iC]->GetPad((ieta-4)+1);
	}
	}*/

      TGraphErrors *grStoch[nPu-1][nSR];
      TGraphErrors *grConst[nPu-1][nSR];
      TGraphErrors *grNoise[nPu-1][nSR];
      for (unsigned iSR(0); iSR<nSR;++iSR){
	
	for (unsigned ipu(0); ipu<( (nPu>1)?(nPu-1):nPu ); ++ipu){//loop on pu
	  double layval[nLayers];
	  double layvalerr[nLayers];
	  for (unsigned iL(0); iL<nLayers;++iL) {
	    layval[iL] = list[iL]; 
	    layvalerr[iL] = 0;
	  }
	  //layval[10-1] = -1;
	  grStoch[ipu][iSR] = new TGraphErrors(nLayers,layval,sigmaStoch[ipu+1][iSR],layvalerr,sigmaStochErr[ipu+1][iSR]);
	  grConst[ipu][iSR] = new TGraphErrors(nLayers,layval,sigmaConst[ipu+1][iSR],layvalerr,sigmaConstErr[ipu+1][iSR]);
	  grNoise[ipu][iSR] = new TGraphErrors(nLayers,layval,sigmaNoise[ipu+1][iSR],layvalerr,sigmaNoiseErr[ipu+1][iSR]);

	  TGraphErrors *gr=0;
	  for (unsigned iP(0);iP<3;++iP){
	    gr = (iP==0) ? grStoch[ipu][iSR] : (iP==1) ? grConst[ipu][iSR] : grNoise[ipu][iSR];
	    if (!gr) continue;
	    if (iP==0) mycR->cd();
	    else if (iP==1) mycC->cd();
	    else mycN->cd();
	    gPad->SetGridy(1);
	    gr->SetLineColor(iSR>=4?iSR+2:iSR+1);
	    gr->SetMarkerColor(iSR>=4?iSR+2:iSR+1);
	    gr->SetMarkerStyle(ipu+21);
	    gr->SetMinimum(0);
	    if (iP<2) gr->SetMaximum(iP==0?1.1:0.13);
	    else gr->SetMaximum(5);
	    if (iP==0) gr->SetTitle(";Layer removed;sampling term (GeV^{#frac{1}{2}})");
	    else if (iP==1) gr->SetTitle(";Layer removed;constant term");
	    else gr->SetTitle(";Layer removed;noise term (GeV)");
	    gr->Draw( (ipu==0 && iSR==0) ? "APEL" : "PEL");
	    label.str("");
	    if (iP==0){
	      label.str("");
	      label << "SR " << iSR;//pu[ipu+1];
	      leg->AddEntry(gr,label.str().c_str(),"P");
	    }
	    //if (ipu==nPu-2){
	    //label.str("");
	      //label << "SR" << iSR;
	      //lat.SetTextSize(0.06);
	      //lat.DrawLatexNDC(0.2,0.85,label.str().c_str());
	    // leg->Draw("same");
	    //}
	    //if (ipu==0) lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");
	  }
	}
      
      }//loop on sr
      mycR->cd();
      leg->Draw("same");
      lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");
      mycR->Update();
      std::ostringstream lsave;
      lsave << plotDir << "/SamplingTerm_vsLayer" << pSuffix;
      if (dovsE) lsave << "_vsE";
      //lsave << "_sr" << iSR;
      mycR->Print((lsave.str()+".pdf").c_str());
      mycR->Print((lsave.str()+".png").c_str());
      
      mycC->cd();
      leg->Draw("same");
      lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");
      mycC->Update();
      lsave.str("");
      lsave << plotDir << "/ConstantTerm_vsLayer" << pSuffix;
      if (dovsE) lsave << "_vsE";
      //lsave << "_sr" << iSR;
      mycC->Print((lsave.str()+".pdf").c_str());
      mycC->Print((lsave.str()+".png").c_str());
      
      mycN->cd();
      leg->Draw("same");
      lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");
      mycN->Update();
      lsave.str("");
      lsave << plotDir << "/NoiseTerm_vsLayer" << pSuffix;
      if (dovsE) lsave << "_vsE";
      //lsave << "_sr" << iSR;
      mycN->Print((lsave.str()+".pdf").c_str());
      mycN->Print((lsave.str()+".png").c_str());
      
      }//loop on scenarios
      
    }//loop on versions
    
    fcalib->Write();

    std::cout << " --min value " << minval << "found for layer " << minlayer << ", using " << list[ic] << std::endl;

    lToRemove.push_back(list[ic]);//minlayer);

  }//loop on Remove vals


  return 0;
  
}//main
