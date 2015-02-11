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

#include "TDRStyle.h"


double zpos(const unsigned l){
  if (l==0) return 3173.9;
  if (l==1) return 3183.65;
  if (l==2) return 3190.95;
  if (l==3) return 3200.4;
  if (l==4) return 3207.7;
  if (l==5) return 3217.15;
  if (l==6) return 3224.45;
  if (l==7) return 3233.9;
  if (l==8) return 3241.2;
  if (l==9) return 3250.65;
  if (l==10) return 3257.95;
  if (l==11) return 3268.45;
  if (l==12) return 3276.85;
  if (l==13) return 3287.35;
  if (l==14) return 3295.75;
  if (l==15) return 3306.25;
  if (l==16) return 3314.65;
  if (l==17) return 3325.15;
  if (l==18) return 3333.55;
  if (l==19) return 3344.05;
  if (l==20) return 3352.45;
  if (l==21) return 3364.35;
  if (l==22) return 3375.05;
  if (l==23) return 3386.95;
  if (l==24) return 3397.65;
  if (l==25) return 3409.55;
  if (l==26) return 3420.25;
  if (l==27) return 3432.15;
  if (l==28) return 3442.85;
  if (l==29) return 3454.75;
  return 0;
};

int logWeightingScan(){//main
  SetTdrStyle();

  //TString plotDir = "/afs/cern.ch/work/a/amagnan/PFCalEEAna/PLOTS/gitV00-02-12/version12/gamma/200um/";
  TString plotDir = "/afs/cern.ch/user/a/amagnan/SLHC/PFCal/PFCalEE/analysis/PLOTS/gitV00-02-12/version12/gamma/200um/";
  std::string pteta = "eta19_et60";

  const unsigned nPu = 1;//2;
  unsigned pu[nPu] = {0};//,140};
  //const double theta = 0.297;//0.361;//0.297;

  const unsigned nScans = 50;
  const double wStart = 1.;
  const double wStep = (6.-wStart)/nScans;
  const unsigned nLayers = 30;

  TFile *fin[nPu];

  TFile *fout = TFile::Open(("PLOTS/LogWeightingStudy_"+pteta+".root").c_str(),"RECREATE");
  fout->mkdir("scan");
  fout->mkdir("scan/xpos");
  fout->mkdir("scan/ypos");
  for (unsigned iS(0); iS<nScans;++iS){
    std::ostringstream lName;
    lName << "scan/xpos/scan_" << wStart+iS*wStep;
    fout->mkdir(lName.str().c_str());
    lName.str("");
    lName << "scan/ypos/scan_" << wStart+iS*wStep;
    fout->mkdir(lName.str().c_str());
  }

  fout->cd();
  TH1F *p_xt = new TH1F("p_xt",";x truth (mm)",100,-10,10); 
  TH1F *p_yt = new TH1F("p_yt",";y truth (mm)",100,-10,10);//200,1170,1370); 


  TH1F *p_posx[nPu][nLayers][nScans];
  TH1F *p_posy[nPu][nLayers][nScans];
  TH1F *p_wx[nPu][nLayers][3];
  TH1F *p_wy[nPu][nLayers][3];

  TH2F *p_Exy[nPu][nLayers];
  TProfile *p_deltavsreco_x[nPu][nLayers];
  TProfile *p_deltavsreco_y[nPu][nLayers];
  TProfile *p_deltavsreco_y_wmin[nPu][nLayers];
  TH2F *p_recovstruth_x[nPu][nLayers];
  TH2F *p_recovstruth_y[nPu][nLayers];


  const unsigned nCanvas = 5;  
  TCanvas *mycx[nCanvas];
  TCanvas *mycy[nCanvas];
  for (unsigned iC(0);iC<nCanvas;++iC){
    std::ostringstream lName;
    lName << "mycx" << iC;
    mycx[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1500,1000);
    mycx[iC]->Divide(3,2);
    lName.str("");
    lName << "mycy" << iC;
    mycy[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1500,1000);
    mycy[iC]->Divide(3,2);
  }

  TCanvas *mycW = new TCanvas("mycW","mycW",1);

  TGraphErrors *grX[nPu][nLayers];
  TGraphErrors *grY[nPu][nLayers];

  double resxmin[nPu][nLayers];
  double wxmin[nPu][nLayers];
  double resymin[nPu][nLayers];
  double wymin[nPu][nLayers];
  double lay[nLayers];

  for (unsigned ipu(0); ipu<nPu; ++ipu){//loop on pu
    std::ostringstream linputStr;
    linputStr << plotDir << "/" << pteta << "_pu" << pu[ipu] ;
    //linputStr << "_logweight";
    linputStr << ".root";
    fin[ipu] = TFile::Open(linputStr.str().c_str());
    if (!fin[ipu]) {
      std::cout << " -- Error, input file " << linputStr.str() << " cannot be opened. Skipping..." << std::endl;
      continue;
    }
    else std::cout << " -- File " << linputStr.str() << " successfully opened." << std::endl;

    TTree *tree = (TTree*)gDirectory->Get("EcellsSR2");
    if (!tree) {
      std::cout << " Tree not found! " << std::endl;
      return 1;
    }
    std::vector<std::vector<double> > Exy;
    std::vector<double> init;
    init.resize(25,0);
    Exy.resize(nLayers,init);
    std::vector<double> truthPosX;
    truthPosX.resize(nLayers,0);
    std::vector<double> truthPosY;
    truthPosY.resize(nLayers,0);

    std::ostringstream label;
    label << "pu" << pu[ipu];
    fout->mkdir(label.str().c_str());
    fout->cd(label.str().c_str());

    for (unsigned iL(0);iL<nLayers;++iL){
      lay[iL] = iL;
      resxmin[ipu][iL] = 100;
      wxmin[ipu][iL] = 10;
      resymin[ipu][iL] = 100;
      wymin[ipu][iL] = 10;
      label.str("");   
      label << "pu" << pu[ipu];
      fout->cd(label.str().c_str());
      label.str("");   
      label << "grX_pu" << pu[ipu] << "_" << iL;
      grX[ipu][iL] = new TGraphErrors();
      grX[ipu][iL]->SetName(label.str().c_str());
      label.str("");   
      label << "grY_pu" << pu[ipu] << "_" << iL;
      grY[ipu][iL] = new TGraphErrors();
      grY[ipu][iL]->SetName(label.str().c_str());

      label.str("");   
      label << "Exy_pu"<< pu[ipu] << "_" << iL;
      p_Exy[ipu][iL] = new TH2F(label.str().c_str(),";x idx;y idx; E (mips)",
				5,0,5,5,0,5);

      label.str("");   
      label << "deltavsreco_x_pu"<< pu[ipu] << "_" << iL;
      p_deltavsreco_x[ipu][iL] = new TProfile(label.str().c_str(),";x reco (mm);x_{reco}-x_{truth} (mm);",
					      30,-15,15,-100,100);
      label.str("");
      label << "deltavsreco_y_pu"<< pu[ipu] << "_" << iL;
      p_deltavsreco_y[ipu][iL] = new TProfile(label.str().c_str(),";y reco (mm);y_{reco}-y_{truth} (mm);",
					      30,-15,15,-100,100);
       label.str("");
      label << "deltavsreco_y_wmin_pu"<< pu[ipu] << "_" << iL;
      p_deltavsreco_y_wmin[ipu][iL] = new TProfile(label.str().c_str(),";y reco (mm);y_{reco}-y_{truth} (mm);",
						   30,-15,15,-100,100);
      label.str("");   
      label << "recovstruth_x_pu"<< pu[ipu] << "_" << iL;
      p_recovstruth_x[ipu][iL] = new TH2F(label.str().c_str(),";x_{truth} (mm);x reco (mm)",
					  30,-15,15,30,-15,15);
      label.str("");
      label << "recovstruth_y_pu"<< pu[ipu] << "_" << iL;
      p_recovstruth_y[ipu][iL] = new TH2F(label.str().c_str(),";y_{truth} (mm);y reco (mm)",
					  30,-15,15,//200,1170,1370,
					  30,-15,15
					  ); 



      label.str("");     
      label << "TruthPosX_" << iL;
      tree->SetBranchAddress(label.str().c_str(),&truthPosX[iL]);
      label.str("");     
      label << "TruthPosY_" << iL;
      tree->SetBranchAddress(label.str().c_str(),&truthPosY[iL]);
      
      for (unsigned i(0);i<5;++i){
	label.str("");     
	label << "wx_pu" << pu[ipu] << "_" << iL << "_" << i;
	p_wx[ipu][iL][i] = new TH1F(label.str().c_str(),
				  ";wx;events",
				  100,-10,0);
	label.str("");     
	label << "wy_pu" << pu[ipu] << "_" << iL << "_" << i;
	p_wy[ipu][iL][i] = new TH1F(label.str().c_str(),
				  ";wy;events",
				  100,-10,0);
      }
      for (unsigned iS(0); iS<nScans;++iS){
	label.str("");
	label << "scan/xpos/scan_" << wStart+iS*wStep;
	fout->cd(label.str().c_str());
	label.str("");
	label << "posx_pu" << pu[ipu] << "_" << iL << "_" << iS;
	p_posx[ipu][iL][iS] = new TH1F(label.str().c_str(),
				     ";x-x_{truth} (mm);events",
				     100,-5,5);
	label.str("");
	label << "scan/ypos/scan_" << wStart+iS*wStep;
	fout->cd(label.str().c_str());
	label.str("");
	label << "posy_pu" << pu[ipu] << "_" << iL << "_" << iS;
	p_posy[ipu][iL][iS] = new TH1F(label.str().c_str(),
				     ";y-y_{truth} (mm);events",
				     100,-5,5);
      }
      
      for (unsigned iy(0);iy<5;++iy){
	for (unsigned ix(0);ix<5;++ix){
	  unsigned idx = 5*iy+ix;
	  label.str("");     
	  label << "E_" << iL << "_" << idx;
	  tree->SetBranchAddress(label.str().c_str(),&Exy[iL][idx]);
	}
      }
    }//loop on layers

    unsigned nEvts = tree->GetEntries();
    for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries

      if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;

      tree->GetEntry(ievt);

      for (unsigned iL(0);iL<nLayers;++iL){
	double Etot = 0;
	double Ex[5] = {0,0,0,0,0};
	double Ey[5] = {0,0,0,0,0};
	for (unsigned idx(0);idx<25;++idx){
	  if (iL>22) Etot += Exy[iL][idx];
	  else if ((idx>5 && idx<9)||
		   (idx>10 && idx<14)||
		   (idx>15 && idx<19)) Etot += Exy[iL][idx];
	}
	if (iL>22){
	  Ex[0] = Exy[iL][0]+Exy[iL][5]+Exy[iL][10]+Exy[iL][15]+Exy[iL][20];
	  Ex[1] = Exy[iL][1]+Exy[iL][6]+Exy[iL][11]+Exy[iL][16]+Exy[iL][21];
	  Ex[2] = Exy[iL][2]+Exy[iL][7]+Exy[iL][12]+Exy[iL][17]+Exy[iL][22];
	  Ex[3] = Exy[iL][3]+Exy[iL][8]+Exy[iL][13]+Exy[iL][18]+Exy[iL][23];
	  Ex[4] = Exy[iL][4]+Exy[iL][9]+Exy[iL][14]+Exy[iL][19]+Exy[iL][24];
	  Ey[0] = Exy[iL][0]+Exy[iL][1]+Exy[iL][2]+Exy[iL][3]+Exy[iL][4];
	  Ey[1] = Exy[iL][5]+Exy[iL][6]+Exy[iL][7]+Exy[iL][8]+Exy[iL][9];
	  Ey[2] = Exy[iL][10]+Exy[iL][11]+Exy[iL][12]+Exy[iL][13]+Exy[iL][14];
	  Ey[3] = Exy[iL][15]+Exy[iL][16]+Exy[iL][17]+Exy[iL][18]+Exy[iL][19];
	  Ey[4] = Exy[iL][20]+Exy[iL][21]+Exy[iL][22]+Exy[iL][23]+Exy[iL][24];
	}
	else {
	  Ex[0] = Exy[iL][6]+Exy[iL][11]+Exy[iL][16];
	  Ex[1] = Exy[iL][7]+Exy[iL][12]+Exy[iL][17];
	  Ex[2] = Exy[iL][8]+Exy[iL][13]+Exy[iL][18];
	  Ey[0] = Exy[iL][6]+Exy[iL][7]+Exy[iL][8];
	  Ey[1] = Exy[iL][11]+Exy[iL][12]+Exy[iL][13];
	  Ey[2] = Exy[iL][16]+Exy[iL][17]+Exy[iL][18];
	}

	double simplex = 0;
	double simpley = 0;
	if (Etot!=0) {
	  if (iL>22) {
	    simplex = 10*(2*Ex[4]+Ex[3]-Ex[1]-2*Ex[0])/Etot;
	    simpley = 10*(2*Ey[4]+Ey[3]-Ey[1]-2*Ey[0])/Etot;
	  }
	  else {
	    simplex = 10*(Ex[2]-Ex[0])/Etot;
	    simpley = 10*(Ey[2]-Ey[0])/Etot;
	  }
	}
	double sign = fabs(truthPosX[iL])/truthPosX[iL];
	unsigned cellCenter = static_cast<unsigned>((fabs(truthPosX[iL])+5)/10.)*10;
	double xt = truthPosX[iL]-sign*cellCenter;
	sign = fabs(truthPosY[iL])/truthPosY[iL];
	cellCenter = static_cast<unsigned>((fabs(truthPosY[iL])+5)/10.)*10;
	double yt=truthPosY[iL]-sign*cellCenter;
	p_xt->Fill(xt);
	p_yt->Fill(yt);
	p_deltavsreco_x[ipu][iL]->Fill(xt,simplex-xt);
	p_deltavsreco_y[ipu][iL]->Fill(yt,simpley-yt);
	p_recovstruth_x[ipu][iL]->Fill(xt,simplex);
	p_recovstruth_y[ipu][iL]->Fill(yt,simpley);

	
	for (unsigned idx(0);idx<25;++idx){
	  p_Exy[ipu][iL]->Fill(idx%5,idx/5,Exy[iL][idx]/Etot);
	}
	double wx[6][nScans];
	double wy[6][nScans];
	  
	for (unsigned i(0);i<6;++i){
	  for (unsigned iS(0); iS<nScans;++iS){
	    wx[i][iS] = 0;
	    wy[i][iS] = 0;
	  }
	}
	for (unsigned i(0);i<5;++i){
	  p_wx[ipu][iL][i]->Fill(log(Ex[i]/Etot));
	  p_wy[ipu][iL][i]->Fill(log(Ey[i]/Etot));
	  for (unsigned iS(0); iS<nScans;++iS){
	    double w0 = wStart+iS*wStep;
	    wx[i][iS] = std::max(0.,log(Ex[i]/Etot)+w0);
	    wy[i][iS] = std::max(0.,log(Ey[i]/Etot)+w0);
	    // if (log(Ex[i]/Etot)+w0<0)
	    //   std::cout << " - iL= " << iL << " i=" << i << " w0=" << w0 
	    // 		<< " logEx=" << log(Ex[i]/Etot)
	    // 		<< " wx " << wx[i][iS] 
	    // 		<< std::endl;
	    // if (log(Ey[i]/Etot)+w0<0) 
	    //   std::cout << " - iL= " << iL << " i=" << i << " w0=" << w0 
	    // 		<< " logEy=" << log(Ey[i]/Etot)
	    // 		<< " wy " << wy[i][iS] 
	    // 		<< std::endl;
	    wx[5][iS] += wx[i][iS];
	    wy[5][iS] += wy[i][iS];
	  }
	}
	for (unsigned iS(0); iS<nScans;++iS){
	  double x = 0;//10*(wx[2][iS]-wx[0][iS])/wx[3][iS];
	  if (wx[5][iS]!=0) {
	    if (iL>22) x = 10*(2*wx[4][iS]+wx[3][iS]-wx[1][iS]-2*wx[0][iS])/wx[5][iS];
	    else x = 10*(wx[2][iS]-wx[0][iS])/wx[5][iS];
	  }
	  double y = 0;//10*(wy[2][iS]-wy[0][iS])/wy[3][iS];
	  if (wy[5][iS]!=0) {
	    if (iL>22) y = 10*(2*wy[4][iS]+wy[3][iS]-wy[1][iS]-2*wy[0][iS])/wy[5][iS];
	    else y = 10*(wy[2][iS]-wy[0][iS])/wy[5][iS];
	  }
	  //if (fabs(y-yt)>5) std::cout << " --- iL=" << iL << " iS=" << iS 
	  //<< " x=" << x << " xt=" << xt 
	  //<< " y=" << y << " yt=" << yt 
	  //<< std::endl;
	  p_posx[ipu][iL][iS]->Fill(x-xt);
	  p_posy[ipu][iL][iS]->Fill(y-yt);
	  if (iS==16) p_deltavsreco_y_wmin[ipu][iL]->Fill(yt,y-yt);

	}

      }//loop on layers
      
    }//loop on entries

    TLatex lat;
    char buf[500];
    for (unsigned iL(0);iL<nLayers;++iL){
      for (unsigned iS(0); iS<nScans;++iS){
	double w0 = wStart+iS*wStep;
	p_posx[ipu][iL][iS]->Fit("gaus","0+","",-1.5,1.5);
	TF1 *fitx = p_posx[ipu][iL][iS]->GetFunction("gaus");
	p_posy[ipu][iL][iS]->Fit("gaus","0+","",-1.5,1.5);
	TF1 *fity = p_posy[ipu][iL][iS]->GetFunction("gaus");
	//grX[ipu][iL]->SetPoint(iS,w0,p_posx[ipu][iL][iS]->GetRMS());
	//grX[ipu][iL]->SetPointError(iS,0,p_posx[ipu][iL][iS]->GetRMSError());
	//grY[ipu][iL]->SetPoint(iS,w0,p_posy[ipu][iL][iS]->GetRMS());
	//grY[ipu][iL]->SetPointError(iS,0,p_posy[ipu][iL][iS]->GetRMSError());
	double xval = 0;
	if (fitx->GetChisquare()/fitx->GetNDF()<20 && fitx->GetParameter(2)<p_posx[ipu][iL][iS]->GetRMS()){
	  xval = fitx->GetParameter(2);
	  grX[ipu][iL]->SetPointError(iS,0,fitx->GetParError(2));
	}
	else {
	  xval = p_posx[ipu][iL][iS]->GetRMS();
	  grX[ipu][iL]->SetPointError(iS,0,p_posx[ipu][iL][iS]->GetRMSError());
	}
	grX[ipu][iL]->SetPoint(iS,w0,xval);
	if (xval < resxmin[ipu][iL]){
	  resxmin[ipu][iL] = xval;
	  wxmin[ipu][iL] = w0;
	}
	double yval = 0;
	if (fity->GetChisquare()/fity->GetNDF()<20 && fity->GetParameter(2)<p_posy[ipu][iL][iS]->GetRMS()){
	  yval = fity->GetParameter(2);
	  grY[ipu][iL]->SetPointError(iS,0,fity->GetParError(2));
	}
	else {
	  yval = p_posy[ipu][iL][iS]->GetRMS();
	  grY[ipu][iL]->SetPointError(iS,0,p_posy[ipu][iL][iS]->GetRMSError());
	}
	grY[ipu][iL]->SetPoint(iS,w0,yval);
	if (yval < resymin[ipu][iL]){
	  resymin[ipu][iL] = yval;
	  wymin[ipu][iL] = w0;
	}
      }
      grX[ipu][iL]->SetTitle(";W0; #sigma(x-xt) (mm)");
      grY[ipu][iL]->SetTitle(";W0; #sigma(y-yt) (mm)");
      grX[ipu][iL]->SetLineColor(ipu+1);
      grX[ipu][iL]->SetMarkerColor(ipu+1);
      grX[ipu][iL]->SetMarkerStyle(ipu+21);
      grY[ipu][iL]->SetLineColor(ipu+1);
      grY[ipu][iL]->SetMarkerColor(ipu+1);
      grY[ipu][iL]->SetMarkerStyle(ipu+21);

      mycx[iL/6]->cd(iL%6+1);
      grX[ipu][iL]->Draw(ipu==0?"APL":"PLsame");
      //gStyle->SetStatX(0.4);
      //gStyle->SetStatY(1.0);
      //p_Exy[ipu][iL]->Draw("colztext");
      sprintf(buf,"Layer %d",iL);
      lat.DrawLatexNDC(0.4,0.85,buf);
      mycy[iL/6]->cd(iL%6+1);
      grY[ipu][iL]->Draw(ipu==0?"APL":"PLsame");
      lat.DrawLatexNDC(0.4,0.85,buf);
 
    }//loop on layers

    //return 1;
  }//loop on pu

  for (unsigned iC(0);iC<nCanvas;++iC){
    mycx[iC]->Update();
    std::ostringstream lsave;
    lsave << "PLOTS/logWeighted_x_" << 6*iC << "_" << 6*iC+5 << "_" << pteta << ".pdf";
    mycx[iC]->Print(lsave.str().c_str());
    mycx[iC]->Update();

    lsave.str("");
    lsave << "PLOTS/logWeighted_y_" << 6*iC << "_" << 6*iC+5 << "_" << pteta << ".pdf";
    mycy[iC]->Print(lsave.str().c_str());
  }

  TGraph *grW[4];
  for (unsigned ipu(0); ipu<nPu; ++ipu){//loop on pu
    std::cout << " --Processing pu " << pu[ipu] << std::endl;

    grW[2*ipu] = new TGraph(nLayers,lay,wxmin[ipu]);
    grW[2*ipu+1] = new TGraph(nLayers,lay,wymin[ipu]);
    for (unsigned iL(0);iL<nLayers;++iL){
      std::cout << " Layer " << iL 
		<< " wmin " << wxmin[ipu][iL] << " " << wymin[ipu][iL] << " " << (wxmin[ipu][iL]+wymin[ipu][iL])/2.
	//		<< " minimum x= " << grX[ipu][iL]->GetYaxis()->GetXmin()
    //		<< " minimum y= " << grY[ipu][iL]->GetYaxis()->GetXmin()
    		<< std::endl;
    }
    mycW->cd();
    grW[2*ipu]->SetTitle(";layer;W0");
    grW[2*ipu]->SetMaximum(6);
    grW[2*ipu]->SetMinimum(0);
    grW[2*ipu]->SetLineColor(2*ipu+1);
    grW[2*ipu]->SetMarkerColor(2*ipu+1);
    grW[2*ipu]->SetMarkerStyle(20+2*ipu+1);
    grW[2*ipu+1]->SetLineColor(2*ipu+2);
    grW[2*ipu+1]->SetMarkerColor(2*ipu+2);
    grW[2*ipu+1]->SetMarkerStyle(20+2*ipu+2);
    if (ipu==0) {
      grW[2*ipu]->Draw("APL");
      grW[2*ipu+1]->Draw("PLsame");
    }
    else {
      grW[2*ipu]->Draw("PLsame");
      grW[2*ipu+1]->Draw("PLsame");
    }
  }
  std::ostringstream lsave;
  lsave << "PLOTS/w0minvsLayers_" << pteta << ".pdf";
  mycW->Update();
  mycW->Print(lsave.str().c_str());

  fout->Write();

  return 0;
}//main
