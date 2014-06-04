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

TPad* plot_ratio(TCanvas *canv, bool up){
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

int plotE(){//main

  const unsigned nSmear = 1;
  TString smearFact[nSmear] = {"No smearing"};//,"1% smearing","2% smearing","3% smearing","4% smearing","5% smearing","7% smearing","10% smearing","15% smearing","20% smearing"};

  //const unsigned nS = 7;
  //TString scenario[nS] = {"0","1","2","3","4","5","6"};
  const unsigned nS = 1;
  std::string scenario[nS] = {
    //"pi-/twiceSampling/GeVCal/MipThresh_0p5/ECALloss/"
    //"pi-/twiceSampling/GeVCal/EarlyDecay/MipThresh_0p5/"
    //"e-/twiceSampling/MipThresh_0p5/"
    "e-/"
  };

  TString pSuffix = "_200um";

  const bool addNoiseTerm = false;

  const unsigned nV = 1;
  TString version[nV] = {"8"};//,"0"};
  
  const unsigned nLayers = 31;//9;//33; //54;

  const bool doFrac = false;
  const bool doShower = false;

  TString pDetector = "Total";//"_FHCAL";

  bool doMIPconv = false;

  double MIPtoGeV = 245.71;//0.914;//41.69;//43.97;//0.92;// 41.98;
  double offset = 431;//-1.04;//-4.3;//-38;//-1.06;
  
  char unitStr[10] = "MIPs";

  const unsigned MAX = 8;
  TString type[MAX];
  Float_t sigmaStoch[nSmear][nS][MAX];
  Float_t sigmaStochErr[nSmear][nS][MAX];
  Float_t sigmaConst[nSmear][nS][MAX];
  Float_t sigmaConstErr[nSmear][nS][MAX];
  Float_t sigmaNoise[nSmear][nS][MAX];
  Float_t sigmaNoiseErr[nSmear][nS][MAX];
  
  std::ostringstream saveName;
  bool isPU = false;
  
  
  unsigned genEn[]={10,15,20,25,30,40,50,60,80,100,150,200,300};
  //60,80};//,100,200,300,
  //500};//,1000,2000};
  //unsigned genEn[]={10,20,30,40,60,80};
  //unsigned genEn[]={40,50,60,80,100,200,300,400,500,1000,2000};
  //unsigned genEn[]={5,10,20,25,50,75,100,125,150,175,200,300,500};
  //unsigned genEn[]={10,25,40,50,60,80};
  //unsigned genEn[]={10,40};
  const unsigned nGenEn=sizeof(genEn)/sizeof(unsigned);
  unsigned rebin[20] = {4,4,4,6,6,
			6,6,8,8,20,
			20,100,100,100,100,
			6,6,6,6,6};
  //unsigned rebin[6] = {12,10,6,6,6,6};

  //canvas so they are created only once
  TCanvas *mycL = 0;
  TCanvas *mycF = 0;
  if (doFrac) {
    mycF = new TCanvas("mycF","mycF",1500,750);
    mycL = new TCanvas("mycL","mycL",1500,1000);
  }
  TCanvas *mycE = new TCanvas("mycE","mycE",1500,1000);
  TCanvas *mycPU = 0;
  if (isPU) mycPU = new TCanvas("mycPU","mycPU",1500,1000);

  if (nGenEn>12)
    mycE->Divide(5,3);
  else if (nGenEn > 10)
    mycE->Divide(4,3);
  else if (nGenEn > 6)
    mycE->Divide(5,2);
  else if (nGenEn > 4)
    mycE->Divide(3,2);
  else if (nGenEn > 2)
    mycE->Divide(2,2);
  else 
    mycE->Divide(nGenEn,1);

  const unsigned nCanvas = 4;  
  TCanvas *myc[nCanvas];
  for (unsigned iC(0);iC<nCanvas;++iC){
    std::ostringstream lName;
    lName << "myc" << iC;
    myc[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1);
  }
  
  TPad *upper = plot_ratio(myc[1], true);
  TPad *lower = plot_ratio(myc[1], false);
  if (!upper || !lower){
    std::cout << " Pb..." << upper << " " << lower << std::endl;
    return 1;
  }

  for (unsigned iV(0); iV<nV;++iV){//loop on versions
    for (unsigned iS(0); iS<nS;++iS){//loop on scenarios
      
      if (scenario[iS].find("PU") != scenario[iS].npos) isPU = true;
      
      TString plotDir = "../PLOTS/version"+version[iV]+"/"+scenario[iS]+"/";
      //plotDir += "noWeights/";
      //TString plotDir = "../PLOTS/version_"+version[iV]+"/scenario_"+scenario[iS]+"/";
      //TString plotDir = "../PLOTS/version_"+version[iV]+"/";
      
      for (unsigned iSm(0); iSm<nSmear; ++iSm){//loop on smear
	unsigned smearOption = iSm;
	
	

	TH1F *p_Efrac[nGenEn][nLayers];
	TH1F *p_Etotal[nGenEn];
	TH1F *p_Ereco[nGenEn];
	TH1F *p_EtotalPU = 0;
	TH1F *p_ErecoPU = 0;
	
	TH2F *p_meanFrac = new TH2F("p_meanFrac",";layer;gen E (GeV);<E_{layer}>/E_{tot}",nLayers,0,nLayers,nGenEn,0,nGenEn);
	TH2F *p_rmsFrac = new TH2F("p_rmsFrac",";layer;gen E (GeV);#sigma(E_{layer}/E_{tot})",nLayers,0,nLayers,nGenEn,0,nGenEn);
	
	
	bool isG4File = false;
	for (unsigned iE(0); iE<nGenEn; ++iE){
	  
	  if (scenario[iS].find("PedroPU")!=scenario[iS].npos) genEn[iE]=iE;
	  
	  std::cout << "- Processing energy : " << genEn[iE] 
		    << std::endl;

	  TFile *inputFile = 0;
	  std::ostringstream linputStr;
	  if (doShower) linputStr << plotDir << "CalibHcalHistos_" << pSuffix << ".root";
	  else linputStr << plotDir << "validation_e" << genEn[iE] << pSuffix << ".root";
	  inputFile = TFile::Open(linputStr.str().c_str());
	  if (!inputFile) {
	    std::cout << " -- Error, input file " << linputStr.str() << " cannot be opened. Exiting..." << std::endl;
	    return 1;
	  }
	  else std::cout << " -- File " << inputFile->GetName() << " sucessfully opened." << std::endl;
	  
	  if (isPU) rebin[iE] = 1;//rebin[iE]*4;
	  
	  TString eStr;
	  eStr += genEn[iE];
	  p_meanFrac->GetYaxis()->SetBinLabel(iE+1,eStr);
	  p_rmsFrac->GetYaxis()->SetBinLabel(iE+1,eStr);
	  
	  if (genEn[iE]==30 && doFrac) {
	    mycL->Divide(6,6);
	    gStyle->SetOptStat(0);
	  }
	  
	  std::ostringstream lName;
	  if (doFrac){
	    for (unsigned iL(0); iL<nLayers; ++iL){
	      //std::cout << " -- Processing layer " << iL << std::endl;
	      lName.str("");
	      lName << "p_EfracReco_" << iL;
	      p_Efrac[iE][iL] = (TH1F*)gDirectory->Get(lName.str().c_str());
	      if (!p_Efrac[iE][iL]) {
		std::cout << " -- ERROR, pointer for histogram Efrac is null for energy " << genEn[iE] << " layer: " << iL << ". Exiting..." << std::endl;
		return 1;
	      }
	      
	      //std::cout << " ---- mean frac = " << p_Efrac[iE][iL]->GetMean() << std::endl;
	      
	      p_meanFrac->SetBinContent(iL+1,iE+1,p_Efrac[iE][iL]->GetMean());
	      p_rmsFrac->SetBinContent(iL+1,iE+1,p_Efrac[iE][iL]->GetRMS());
	      
	      if (genEn[iE]==30 && iL<36) {
		mycL->cd(iL+1);
		char buf[500];
		sprintf(buf,"E=%d GeV, layer %d",genEn[iE],iL);
		p_Efrac[iE][iL]->SetTitle(buf);
		p_Efrac[iE][iL]->Draw();
	      }
	      
	    }//loop on layers

	    
	    if (genEn[iE]==30){
	      saveName.str("");
	      saveName << plotDir << "/SimEfraction_" << genEn[iE] << "GeV" ;
	      mycL->Update();
	      mycL->Print((saveName.str()+".png").c_str());
	      mycL->Print((saveName.str()+".pdf").c_str());
	    }
	  }//doFrac

	  lName.str("");
	  if (doShower) lName << "p_EshowerCor_" << genEn[iE];
	  //if (doShower) lName << "p_Eshower_" << genEn[iE] << "_" << 5;
	  else lName << "p_Esim" << pDetector;
	  p_Etotal[iE] = (TH1F*)gDirectory->Get(lName.str().c_str());
	  if (!p_Etotal[iE]){
	    std::cout << " -- ERROR, pointer for histogram " << lName.str() << " is null. Exiting..." << std::endl;
	    return 1;
	  }
	  p_Etotal[iE]->Sumw2();
	  //remove e/pi to recalculate proper one....

	  std::cout << " --- Sim E = entries " << p_Etotal[iE]->GetEntries() 
		    << " mean " << p_Etotal[iE]->GetMean() 
		    << " rms " << p_Etotal[iE]->GetRMS() 
		    << " overflows " << p_Etotal[iE]->GetBinContent(p_Etotal[iE]->GetNbinsX()+1)
		    << std::endl;

	  p_Etotal[iE]->Rebin(rebin[iE]);
	  //p_Etotal[iE]->Rebin(4);

	  lName.str("");
	  lName << "p_Ereco_" << genEn[iE] << "_smear" << smearOption;
	  p_Ereco[iE] = (TH1F*)gDirectory->Get(lName.str().c_str());
	  if (!p_Ereco[iE]){
	    std::cout << " -- ERROR, pointer for histogram Ereco is null. Running on G4 file before digitizer." << std::endl;
	    isG4File = true;
	  }
	  else {
	    std::cout << " --- Reco E = entries " << p_Ereco[iE]->GetEntries() 
		      << " mean " << p_Ereco[iE]->GetMean() 
		      << " rms " << p_Ereco[iE]->GetRMS() 
		      << " overflows " << p_Ereco[iE]->GetBinContent(p_Etotal[iE]->GetNbinsX()+1)
		      << std::endl;

	    p_Ereco[iE]->Rebin(rebin[iE]);
	    //p_Ereco[iE]->Rebin(rebin);
	  }

	}//loop on energies

	//patch for HCAL only setup
	if ((version[iV] == "21" || version[iV] == "22") && pDetector == "ECAL") pDetector = "HCAL";
	if (version[iV] == "21" && pDetector == "Total") pDetector = "SiSci";
	if (version[iV] == "23" && pDetector == "Total") pDetector = "Calice";
	pDetector += pSuffix;

	if (isPU){//isPU
	  p_EtotalPU = new TH1F("p_EtotalPU",";Etot (MeV)",100,p_Etotal[0]->GetBinLowEdge(1),p_Etotal[nGenEn-1]->GetBinLowEdge(p_Etotal[nGenEn-1]->GetNbinsX()+1));
	  p_EtotalPU->Sumw2();
	  for (unsigned iE(0); iE<nGenEn; ++iE){
	    for (int iB(1); iB<p_Etotal[iE]->GetNbinsX()+1;++iB){//loop on bins
	      p_EtotalPU->Fill(p_Etotal[iE]->GetBinCenter(iB),p_Etotal[iE]->GetBinContent(iB));
	    }
	    if (!isG4File){
	      if (iE==0) {
		p_ErecoPU = new TH1F("p_ErecoPU",";Etot (MIPs)",100,p_Ereco[0]->GetBinLowEdge(1),p_Ereco[nGenEn-1]->GetBinLowEdge(p_Ereco[nGenEn-1]->GetNbinsX()+1));
		p_ErecoPU->Sumw2();
	      }
	      for (int iB(1); iB<p_Ereco[iE]->GetNbinsX()+1;++iB){//loop on bins
		p_ErecoPU->Fill(p_Ereco[iE]->GetBinCenter(iB),p_Ereco[iE]->GetBinContent(iB));
	      }
	    }
	  }
  
	  if (!isG4File){
	    mycPU->Clear();
	    mycPU->cd();
	    gStyle->SetOptStat(1111110);
	    p_ErecoPU->Draw();
	  
	    saveName.str("");
	    saveName << plotDir << "/PUTotalRecoE_smear" << smearOption;
	    mycPU->Update();
	    mycPU->Print((saveName.str()+".png").c_str());
	    mycPU->Print((saveName.str()+".pdf").c_str());
	  
	  }
	  mycPU->Clear();
	  mycPU->cd();
	  gStyle->SetOptStat(1111110);
	  p_EtotalPU->Draw();
	
	  saveName.str("");
	  saveName << plotDir << "/PUTotalSimE" ;
	  mycPU->Update();
	  mycPU->Print((saveName.str()+".png").c_str());
	  mycPU->Print((saveName.str()+".pdf").c_str());
	}//isPU

	if (doFrac){
	  //draw energy fractions
	  gStyle->SetOptStat(0);
	  mycF->Clear();
	  mycF->Divide(2,1);
	  mycF->cd(1);
	  p_meanFrac->Draw("colz");
	  mycF->cd(2);
	  p_rmsFrac->Draw("colz");
      
	  saveName.str("");
	  saveName << plotDir << "/SimEfractionIntegrated";
	  mycF->Update();
	  mycF->Print((saveName.str()+".png").c_str());
	  mycF->Print((saveName.str()+".pdf").c_str());
	}
      
	gStyle->SetOptStat(0);

	//draw calibration curves
	TGraphErrors *calib = new TGraphErrors();
	calib->SetName("calib");
	calib->SetMarkerStyle(21);
	calib->SetTitle("");
	TGraphErrors *reso = (TGraphErrors *) calib->Clone("reso");
	TGraphErrors *calibFit = (TGraphErrors *) calib->Clone("calibFit");
	TGraphErrors *deltaFit = (TGraphErrors *) calib->Clone("deltaFit");
	TGraphErrors *resoFit = (TGraphErrors *) calib->Clone("resoFit");
	TGraphErrors *calibReco = (TGraphErrors *) calib->Clone("calibReco");
	calibReco->SetMarkerStyle(22);
	calibReco->SetMarkerColor(6);
	calibReco->SetLineColor(6);
	TGraphErrors *resoReco = (TGraphErrors *) calibReco->Clone("resoReco");
	TGraphErrors *calibRecoFit = (TGraphErrors *) calibReco->Clone("calibRecoFit");
	TGraphErrors *resoRecoFit = (TGraphErrors *) calibReco->Clone("resoRecoFit");
  
	//simhits
	for (unsigned iE(0); iE<nGenEn; ++iE){
	  std::cout << "- Processing energy : " << genEn[iE] << std::endl;
	  //plot total E
	  mycE->cd(iE+1);
	  gStyle->SetOptFit(0);
	  p_Etotal[iE]->GetXaxis()->SetRangeUser(p_Etotal[iE]->GetMean()-10*p_Etotal[iE]->GetRMS(),p_Etotal[iE]->GetMean()+10*p_Etotal[iE]->GetRMS());
	  p_Etotal[iE]->Draw();
	  if (pDetector == "HCAL" || pDetector == "Total" || pDetector == "SiSci") 
	    p_Etotal[iE]->Fit("gaus","LR+","",
			      p_Etotal[iE]->GetMean()-2*p_Etotal[iE]->GetRMS(),
			      p_Etotal[iE]->GetMean()+2*p_Etotal[iE]->GetRMS());

	  else p_Etotal[iE]->Fit("gaus","LR+","",
				 p_Etotal[iE]->GetMean()-2*p_Etotal[iE]->GetRMS(),
				 p_Etotal[iE]->GetMean()+2*p_Etotal[iE]->GetRMS());


	  TF1 *fitResult = p_Etotal[iE]->GetFunction("gaus");
	  TLatex lat;
	  double latx = std::max(0.,p_Etotal[iE]->GetMean()-5*p_Etotal[iE]->GetRMS());
	  double laty = p_Etotal[iE]->GetMaximum();
	  char buf[500];
	  sprintf(buf,"<E> = %3.3f %s",p_Etotal[iE]->GetMean(),doMIPconv?"MIP":unitStr);
	  lat.DrawLatex(latx,laty*0.9,buf);
	  sprintf(buf,"RMS = %3.3f #pm %3.1f %s",p_Etotal[iE]->GetRMS(),p_Etotal[iE]->GetRMSError(),doMIPconv?"MIP":unitStr);
	  lat.DrawLatex(latx,laty*0.8,buf);
	  sprintf(buf,"RMS/mean = %3.3f",p_Etotal[iE]->GetRMS()/p_Etotal[iE]->GetMean());
	  lat.DrawLatex(latx,laty*0.7,buf);
	  sprintf(buf,"<Efit> = %3.3f +/- %3.3f %s",fitResult->GetParameter(1),fitResult->GetParError(1),doMIPconv?"MIP":unitStr);
	  lat.DrawLatex(latx,laty*0.6,buf);
	  sprintf(buf,"RMSfit = %3.3f +/- %3.3f %s",fitResult->GetParameter(2),fitResult->GetParError(2),doMIPconv?"MIP":unitStr);
	  lat.DrawLatex(latx,laty*0.5,buf);
	  sprintf(buf,"RMS/meanfit = %3.3f",fitResult->GetParameter(2)/fitResult->GetParameter(1));
	  lat.DrawLatex(latx,laty*0.4,buf);

	  sprintf(buf,"chi2/NDF = %3.3f/%d = %3.3f",fitResult->GetChisquare(),fitResult->GetNDF(),fitResult->GetChisquare()/fitResult->GetNDF());
          lat.DrawLatex(latx,laty*0.3,buf);


	  Int_t np=calib->GetN();
	  calib->SetPoint(np,genEn[iE],p_Etotal[iE]->GetMean());
	  calib->SetPointError(np,0.0,p_Etotal[iE]->GetMeanError());
	  reso->SetPoint(np,1/sqrt(genEn[iE]),p_Etotal[iE]->GetRMS()/p_Etotal[iE]->GetMean());
	  reso->SetPointError(np,0,p_Etotal[iE]->GetRMSError()/p_Etotal[iE]->GetMean());
	  calibFit->SetPoint(np,genEn[iE],fitResult->GetParameter(1));
	  calibFit->SetPointError(np,0.0,fitResult->GetParError(1));
	  deltaFit->SetPoint(np,genEn[iE],( ((fitResult->GetParameter(1)-offset)/MIPtoGeV)-genEn[iE])/genEn[iE]);
	  deltaFit->SetPointError(np,0.0,fitResult->GetParError(1)/MIPtoGeV*1./genEn[iE]);
	  resoFit->SetPoint(np,1/sqrt(genEn[iE]),fitResult->GetParameter(2)/fitResult->GetParameter(1));
	  resoFit->SetPointError(np,0,fitResult->GetParError(2)/fitResult->GetParameter(1));

	}//loop on energies


	saveName.str("");
	saveName << plotDir << "/SimG4Etotal_" << pDetector;
	if (doShower) saveName << "_Shower";
	mycE->Update();
	mycE->Print((saveName.str()+".png").c_str());
	mycE->Print((saveName.str()+".pdf").c_str());

	//return 1;

	if (!isG4File){
	  //recohits
	  for (unsigned iE(0); iE<nGenEn; ++iE){
	    //plot reco E
	    std::cout << "- Processing energy : " << genEn[iE] << std::endl;
	    mycE->cd(iE+1);
	    gStyle->SetOptFit(0);
	    //p_Ereco[iE]->GetXaxis()->SetRangeUser(p_Ereco[iE]->GetMean()-5*p_Ereco[iE]->GetRMS(),p_Ereco[iE]->GetMean()+5*p_Ereco[iE]->GetRMS());
	    p_Ereco[iE]->Draw();
	    p_Ereco[iE]->Fit("gaus");
	    TF1 *fitResult = p_Ereco[iE]->GetFunction("gaus");
	    TLatex lat;
	    char buf[500];
	    sprintf(buf,"<E> = %3.1f MIPs",p_Ereco[iE]->GetMean());
	    lat.DrawLatex(p_Ereco[iE]->GetMean()-5*p_Ereco[iE]->GetRMS(),p_Ereco[iE]->GetMaximum()*0.9,buf);
	    sprintf(buf,"RMS = %3.1f MIPs",p_Ereco[iE]->GetRMS());
	    lat.DrawLatex(p_Ereco[iE]->GetMean()-5*p_Ereco[iE]->GetRMS(),p_Ereco[iE]->GetMaximum()*0.8,buf);
	    sprintf(buf,"<Efit> = %3.1f MIPs",fitResult->GetParameter(1));
	    lat.DrawLatex(p_Ereco[iE]->GetMean()+p_Ereco[iE]->GetRMS(),p_Ereco[iE]->GetMaximum()*0.9,buf);
	    sprintf(buf,"RMSfit = %3.1f MIPs",fitResult->GetParameter(2));
	    lat.DrawLatex(p_Ereco[iE]->GetMean()+p_Ereco[iE]->GetRMS(),p_Ereco[iE]->GetMaximum()*0.8,buf);
      
      
	    Int_t np=calibReco->GetN();
	    calibReco->SetPoint(np,genEn[iE],p_Ereco[iE]->GetMean());
	    calibReco->SetPointError(np,0.0,p_Ereco[iE]->GetMeanError());
	    resoReco->SetPoint(np,1/sqrt(genEn[iE]),p_Ereco[iE]->GetRMS()/p_Ereco[iE]->GetMean());
	    resoReco->SetPointError(np,0,p_Ereco[iE]->GetRMSError()/p_Ereco[iE]->GetMean());
	    calibRecoFit->SetPoint(np,genEn[iE],fitResult->GetParameter(1));
	    calibRecoFit->SetPointError(np,0.0,fitResult->GetParError(1));
	    resoRecoFit->SetPoint(np,1/sqrt(genEn[iE]),fitResult->GetParameter(2)/fitResult->GetParameter(1));
	    resoRecoFit->SetPointError(np,0,fitResult->GetParError(2)/fitResult->GetParameter(1));
      
	  }//loop on energies

	  saveName.str("");
	  saveName << plotDir << "/DigiEreco_smear" << smearOption << "_" << pDetector;
	  mycE->Update();
	  mycE->Print((saveName.str()+".png").c_str());
	  mycE->Print((saveName.str()+".pdf").c_str());
	}

	//draw calib

	const unsigned imax = isG4File ? 4 : 8;



	for(unsigned i=0; i<imax ; i++)
	  {
	    if (i==1) upper->cd();
	    else myc[i%4]->cd();
	    type[i] = i==0 ? "calib" : i==1 ? "calibFit" : i==2 ? "reso" : "resoFit";
	    if (!isG4File && i>3) type[i] = i==4 ? "calibReco" : i==5 ? "calibRecoFit" : i==6 ? "resoReco" : "resoRecoFit";

	    if (nSmear > 1){
	      type[i] += "_smear";
	      type[i] += smearOption;
	    }
	    type[i] += "_"+pDetector;

	    std::cout << "- Processing type : " << type[i] << std::endl;

	    TGraphErrors * gr =( i==0 ? calib : i==1 ? calibFit : i==2 ? reso : resoFit);

	    gr->GetXaxis()->SetLabelSize(0.06);
	    gr->GetXaxis()->SetTitleSize(0.06);
	    gr->GetYaxis()->SetLabelSize(0.06);
	    gr->GetYaxis()->SetTitleSize(0.06);
	    gr->GetXaxis()->SetTitleOffset(0.7);
	    gr->GetYaxis()->SetTitleOffset(0.8);

	    if (!isG4File && i>3) gr = i==4 ? calibReco : i==5 ? calibRecoFit : i==6 ? resoReco : resoRecoFit;
	    if (i==3 && nSmear>1) gr->SetTitle(smearFact[iSm]);
	    else gr->SetTitle("");
	    gr->Draw(i<4? "ap" : "p");
	    gr->GetYaxis()->SetRangeUser(0,gr->GetYaxis()->GetXmax());
	    
	    if(i<2|| (!isG4File && (i==4 || i==5))) { 
	      if (i!=1) gr->GetXaxis()->SetTitle("Beam energy [GeV]");
	      else gr->GetXaxis()->SetTitle("");
	      if (i>2 || (doMIPconv&&i>0)) gr->GetYaxis()->SetTitle("Average energy deposited [MIPs]");
	      else gr->GetYaxis()->SetTitle("Average energy deposited ["+TString(unitStr)+"]"); 
	    }
	    else { 
	      gr->GetXaxis()->SetTitle("1/#sqrt{Beam energy} [1/#sqrt{GeV}]"); 
	      gr->GetYaxis()->SetTitle("Relative energy resolution");     
	    }
	    char buf[500];
	    if(i<2 || (!isG4File && (i==4 || i==5))) {
	      TF1 *fitFunc=new TF1("calib","[0]+[1]*x",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
	      if (i<4) fitFunc->SetLineColor(1);
	      else fitFunc->SetLineColor(6);
	      gr->Fit(fitFunc,"RME");
	      TLatex lat;
	      if (i>3) lat.SetTextColor(6);
	      else lat.SetTextColor(1);
	      sprintf(buf,"<E> #propto a + b #times E ");
	      if (i<2) lat.DrawLatex(genEn[0],gr->GetYaxis()->GetXmax()*0.9,buf);
	      sprintf(buf,"a = %3.3f #pm %3.3f ",fitFunc->GetParameter(0),fitFunc->GetParError(0));
	      lat.DrawLatex(genEn[0]+i/2*150,gr->GetYaxis()->GetXmax()*(0.8-i/2*0.25),buf);
	      sprintf(buf,"b = %3.3f #pm %3.3f",fitFunc->GetParameter(1),fitFunc->GetParError(1));
	      lat.DrawLatex(genEn[0]+i/2*150,gr->GetYaxis()->GetXmax()*(0.7-i/2*0.25),buf);
	      sprintf(buf,"chi2/NDF = %3.3f/%d = %3.3f",fitFunc->GetChisquare(),fitFunc->GetNDF(),fitFunc->GetChisquare()/fitFunc->GetNDF());
	      lat.DrawLatex(genEn[0]+i/2*150,gr->GetYaxis()->GetXmax()*(0.6-i/2*0.25),buf);

	      //draw deltaE/E vs E
	      if (i==1){
		lower->cd();
		gPad->SetLogx(0);
		gPad->SetGridx(1);
		gPad->SetGridy(1);
		TGraphErrors * grDelta = deltaFit;
		grDelta->SetTitle("");
		grDelta->SetMinimum(-0.1);
		grDelta->SetMaximum(0.1);
		grDelta->GetXaxis()->SetLabelSize(0.15);
		grDelta->GetXaxis()->SetTitleSize(0.15);
		grDelta->GetYaxis()->SetLabelSize(0.12);
		grDelta->GetYaxis()->SetTitleSize(0.15);
		grDelta->GetXaxis()->SetTitleOffset(0.5);
		grDelta->GetYaxis()->SetTitleOffset(0.3);

		grDelta->Draw("ap");
		//grDelta->GetYaxis()->SetRangeUser(0,grDelta->GetYaxis()->GetXmax());
		grDelta->GetXaxis()->SetTitle("Beam energy [GeV]");
		grDelta->GetYaxis()->SetTitle("(#Delta E)/E");
	      }

	    }
	    else
	      {
		//TF1 *fitFunc2 =new TF1("reso","sqrt([0]*x*x+[1]+[2]*x*x*x*x)",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
		TF1 *fitFunc2;
		fitFunc2 =new TF1("reso","sqrt([0]*x*x+[1])",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
		if (addNoiseTerm) fitFunc2 =new TF1("reso","sqrt([0]*x*x+[1]+[2]*x*x*x*x)",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
		//if (i<4) fitFunc2 =new TF1("reso","sqrt([0]*x*x+[1])",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
		//else fitFunc2 =new TF1("reso","sqrt([0]*x*x+[1])",gr->GetXaxis()->GetXmin(),0.2);//gr->GetXaxis()->GetXmax());

		fitFunc2->SetParameter(0,0.2);
		fitFunc2->SetParLimits(0,0,1);
		fitFunc2->SetParameter(1,0);
		fitFunc2->SetParLimits(1,0,1);
		if (addNoiseTerm) {
		  fitFunc2->SetParameter(2,0.18*0.18);
		  fitFunc2->SetParLimits(2,0,0);
		  fitFunc2->FixParameter(2,0.18*0.18);
		}
		if (pDetector == "HCAL" || pDetector == "Total" || pDetector == "SiSci") {
		  fitFunc2->SetParameter(0,0.5);
		  fitFunc2->SetParLimits(0,0,2);
		  fitFunc2->SetParameter(1,0);
		  fitFunc2->SetParLimits(1,0,1);
		  if (addNoiseTerm) {
		    fitFunc2->SetParameter(2,0.18*0.18);
		    fitFunc2->SetParLimits(2,0,0);
		    fitFunc2->FixParameter(2,0.18*0.18);
		  }
		}
		if (i<4) {
		  //fitFunc2->SetParameter(2,0.);
		  //fitFunc2->SetParLimits(2,0,0);
		  fitFunc2->SetLineColor(1);
		}
		else {
		  //fitFunc2->SetLineColor(6);
		  //fitFunc2->SetParameter(2,0.);
		  //fitFunc2->SetParLimits(2,0,0);
		}
		gr->Fit(fitFunc2,"RME");
		sigmaStoch[iSm][iS][i] = sqrt(fitFunc2->GetParameter(0));
		sigmaStochErr[iSm][iS][i] = fitFunc2->GetParError(0)/(2*sigmaStoch[iSm][iS][i]);
		sigmaConst[iSm][iS][i] = sqrt(fitFunc2->GetParameter(1));
		sigmaConstErr[iSm][iS][i] = fitFunc2->GetParError(1)/(2*sigmaConst[iSm][iS][i]);
		if (addNoiseTerm) {
		  sigmaNoise[iSm][iS][i] = sqrt(fitFunc2->GetParameter(2));
		  sigmaNoiseErr[iSm][iS][i] = fitFunc2->GetParError(2)/(2*sigmaNoise[iSm][iS][i]);
		}
		TLatex lat;
		if (i>3) lat.SetTextColor(6);
		else lat.SetTextColor(1);
		if (addNoiseTerm) sprintf(buf,"#frac{#sigma}{E} #propto #frac{s}{#sqrt{E}} #oplus c #oplus #frac{n}{E}");
		else sprintf(buf,"#frac{#sigma}{E} #propto #frac{s}{#sqrt{E}} #oplus c");
		double Emin = 1/sqrt(genEn[nGenEn-1]);
		if (i<4) lat.DrawLatex(Emin,gr->GetYaxis()->GetXmax()*0.9,buf);
		sprintf(buf,"s=%3.3f #pm %3.3f",sigmaStoch[iSm][iS][i],sigmaStochErr[iSm][iS][i]);
		lat.DrawLatex(Emin,gr->GetYaxis()->GetXmax()*(0.9-i/2*0.1-i/4*0.3),buf);
		sprintf(buf,"c=%3.3f #pm %3.3f",sigmaConst[iSm][iS][i],sigmaConstErr[iSm][iS][i]);
		lat.DrawLatex(Emin,gr->GetYaxis()->GetXmax()*(0.8-i/2*0.1-i/4*0.3),buf);
		sprintf(buf,"chi2/NDF = %3.3f/%d = %3.3f",fitFunc2->GetChisquare(),fitFunc2->GetNDF(),fitFunc2->GetChisquare()/fitFunc2->GetNDF());
		lat.DrawLatex(Emin,gr->GetYaxis()->GetXmax()*(0.4-i/2*0.1-i/4*0.3),buf);
		//if (i>3){
		if (addNoiseTerm) {
		  sprintf(buf,"n=%3.3f #pm %3.3f",sigmaNoise[iSm][iS][i],sigmaNoiseErr[iSm][iS][i]);
		  lat.DrawLatex(Emin,gr->GetYaxis()->GetXmax()*(0.7-i/2*0.1-i/4*0.3),buf);
		}
	      }
	    myc[i%4]->Update();
	    if (doShower) {
	      myc[i%4]->Print(plotDir+"/"+type[i]+"_Shower.pdf");
	      myc[i%4]->Print(plotDir+"/"+type[i]+"_Shower.png");
	    }
	    else {
	      myc[i%4]->Print(plotDir+"/"+type[i]+".pdf");
	      myc[i%4]->Print(plotDir+"/"+type[i]+".png");
	    }
	  }

      }//loop on smear options


    }//loop on scenarios

  }//loop on versions
  
  for (unsigned iV(0); iV<nV;++iV){
    std::cout << "version " << version[iV] << std::endl;
    std::cout << "scenario & type & sigmaStoch & sigmaConst";
    if (addNoiseTerm) std::cout << " & sigmaNoise";
    std::cout << " \\\\ \n" 
	      << "\\hline\n";
    for (unsigned iS(0); iS<nS;++iS){
      for (unsigned i(2); i<MAX;++i){
	for (unsigned iSm(0); iSm<nSmear; ++iSm){//loop on smear
	  if (i==4 || i==5) continue;
	  std::cout << scenario[iS] << " & " << type[i] << " & " <<  std::setprecision(3)
		    << "$" << sigmaStoch[iSm][iS][i] << "\\pm" << sigmaStochErr[iSm][iS][i] << "$ & "
		    << "$" << sigmaConst[iSm][iS][i] << "\\pm" << sigmaConstErr[iSm][iS][i] << "$";
	  if (addNoiseTerm) std::cout << " & $" << sigmaNoise[iSm][iS][i] << "\\pm" << sigmaNoiseErr[iSm][iS][i] << "$";
		  std::cout << "\n";
	}
	std::cout <<"\\hline\n";
      }
    }
  }

  if (nSmear > 1){

    mycF->Clear();
    mycF->cd();
    double x[10] = {0,1,2,3,4,5,7,10,15,20};
    double y[nSmear];
    double xerr[nSmear];
    double yerr[nSmear];
    
    double s[nSmear];
    double serr[nSmear];
    
    double lShift = 0.175;

    for (unsigned iSm(0); iSm<nSmear; ++iSm){//loop on smear
      xerr[iSm] = 0;
      y[iSm] = sigmaConst[iSm][0][7];
      yerr[iSm] = sigmaConstErr[iSm][0][7];
      //offset to have on same plot
      s[iSm] = sigmaStoch[iSm][0][7]-lShift;
      serr[iSm] = sigmaStochErr[iSm][0][7];
    }
    TGraphErrors *gr = new TGraphErrors(nSmear,x,y,xerr,yerr);
    TGraphErrors *grs = new TGraphErrors(nSmear,x,s,xerr,serr);
    
    gr->GetXaxis()->SetTitle("Smearing factor (%)");
    gr->GetYaxis()->SetTitle("Constant term");
    gr->SetTitle("Single electrons in HGCAL-EE");
    
    gr->SetMarkerStyle(20);
    
    gr->SetMarkerColor(1);
    gr->SetLineColor(1);
    
    gr->SetMinimum(0);
    gr->SetMaximum(0.04);
    gr->Draw("AP");
    
    TF1 *BE = new TF1("BE","sqrt([0]*[0] + pow(x/100.*1/sqrt([1]),2))",0,20);
    BE->SetParameters(y[0],100);
    BE->SetParLimits(0,1,1);
    BE->SetLineColor(1);
    gr->Fit("BE");
    
    //grBE->Draw("P");
    
    grs->SetMarkerStyle(22);
    grs->SetMarkerColor(2);
    grs->SetLineColor(2);
    grs->Draw("P");
    grs->Fit("pol0");
    TF1 *pol0 = grs->GetFunction("pol0");
    
    TGaxis *ys = new TGaxis(22,0,22,0.04,lShift,lShift+0.04,510,"+L");
    
    ys->SetTextColor(2);
    ys->SetLabelColor(2);
    ys->SetLineColor(2);
    ys->SetTitle("Sampling term");
    ys->Draw();
    
    TLatex lat;
    lat.SetTextColor(2);
    char buf[500];
    sprintf(buf,"<s> = %1.4f #pm %1.4f",pol0->GetParameter(0)+lShift,pol0->GetParError(0));
    lat.DrawLatex(3,0.037,buf);
    
    lat.SetTextColor(1);
    sprintf(buf,"c #propto c_{0} #oplus #frac{x}{#sqrt{n}}, n=%3.0f #pm %3.0f",BE->GetParameter(1),BE->GetParError(1));
    lat.DrawLatex(8,0.005,buf);
    
    //TLegend *leg = new TLegend(0.5,0.12,0.89,0.3);
    //leg->SetFillColor(10);
    //leg->AddEntry(gr,"Standalone simulation","P");
    //leg->Draw("same");
    
    
    saveName.str("");
    saveName << "../PLOTS/version" << version[0] << "/" << scenario[0] << "/Intercalibration";
    mycF->Update();
    mycF->Print((saveName.str()+".png").c_str());
    mycF->Print((saveName.str()+".pdf").c_str());
  }

  return 0;
  
  
}//main
