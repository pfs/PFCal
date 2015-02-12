#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>

#include "TFile.h"
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

int plotMomentumShowerMax(){


  const unsigned nFiles = 1;//30;


  std::ostringstream lstr;
  lstr << "PLOTS/output_layer30_ftfp.root";
  TFile *outputFile = TFile::Open(lstr.str().c_str(),"RECREATE");
  
  if (!outputFile) {
    std::cout << " -- Error, output file " << lstr.str() << " cannot be opened. Please create output directory. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
  }
  outputFile->cd();

  outputFile->mkdir("Si");
  outputFile->mkdir("pipm");
  outputFile->mkdir("protons");
  outputFile->mkdir("neutrons");
  outputFile->cd("protons");
  TH1F *p_logmomentum_p  = new TH1F("p_logmomentum_p",";log(p) (log(GeV));particles",100,-2,2);
  TH1F *p_momentum_p  = new TH1F("p_momentum_p",";p (GeV);particles",1000,0,100);
  TH1F *p_time_p  = new TH1F("p_time_p",";time (ns);particles",1000,0,1000);
  TH2F *p_timevslogp_p = new TH2F("p_timevslogp_p",";log(p) (log(GeV));time (ns);particles",100,-2,2,1000,0,1000);
  TH2F *p_timevsp_p = new TH2F("p_timevsp_p",";p (GeV);time (ns);particles",1000,0,100,1000,0,1000);

  outputFile->cd();
  outputFile->cd("neutrons");
  TH1F *p_logmomentum_n  = new TH1F("p_logmomentum_n",";log(p) (log(GeV));particles",100,-2,2);
  TH1F *p_momentum_n  = new TH1F("p_momentum_n",";p (GeV);particles",1000,0,100);
  TH1F *p_time_n  = new TH1F("p_time_n",";time (ns);particles",1000,0,1000);
  TH2F *p_timevslogp_n = new TH2F("p_timevslogp_n",";log(p) (log(GeV));time (ns);particles",100,-2,2,1000,0,1000);
  TH2F *p_timevsp_n = new TH2F("p_timevsp_n",";p (GeV);time (ns);particles",1000,0,100,1000,0,1000);

  outputFile->cd();
  outputFile->cd("pipm");
  TH1F *p_logmomentum_pi  = new TH1F("p_logmomentum_pi",";log(p) (log(GeV));particles",100,-2,2);
  TH1F *p_momentum_pi  = new TH1F("p_momentum_pi",";p (GeV);particles",1000,0,100);
  TH1F *p_time_pi  = new TH1F("p_time_pi",";time (ns);particles",1000,0,1000);
  TH2F *p_timevslogp_pi = new TH2F("p_timevslogp_pi",";log(p) (log(GeV));time (ns);particles",100,-2,2,1000,0,1000);
  TH2F *p_timevsp_pi = new TH2F("p_timevsp_pi",";p (GeV);time (ns);particles",1000,0,100,1000,0,1000);

  outputFile->cd();
  outputFile->cd("Si");
  TH1F *p_logmomentum_Si  = new TH1F("p_logmomentum_Si",";log(p) (log(GeV));particles",100,-5,2);
  TH1F *p_momentum_Si  = new TH1F("p_momentum_Si",";p (GeV);particles",1000,0,1);
  TH1F *p_time_Si  = new TH1F("p_time_Si",";time (ns);particles",1000,0,1000);
  TH2F *p_timevslogp_Si = new TH2F("p_timevslogp_Si",";log(p) (log(GeV));time (ns);particles",100,-5,2,1000,0,1000);
  TH2F *p_timevsp_Si = new TH2F("p_timevsp_Si",";p (GeV);time (ns);particles",1000,0,1,1000,0,1000);
  TH1F *p_ionisingE = new TH1F("p_ionisingE",";E ionising (MeV);ions",1000,0,1);
  TH2F *p_trkStatus = new TH2F("p_trkStatus",";process name;trk status;ions",10,0,10,6,0,6);
  p_trkStatus->GetXaxis()->SetBinLabel(1,"hadElastic");
  p_trkStatus->GetXaxis()->SetBinLabel(2,"NeutronInelastic");
  p_trkStatus->GetXaxis()->SetBinLabel(3,"ProtonInelastic");
  p_trkStatus->GetXaxis()->SetBinLabel(4,"PionPlusInelastic");

  TH1F *p_stepLength = new TH1F("p_stepLength",";step length (mm);ions",101,0,0.101);
 
  std::map<std::string,unsigned> procnameMap;
  std::map<int,unsigned> pdgIdMap;
  std::pair<std::map<int,unsigned>::iterator,bool> isInserted;
  std::pair<std::map<std::string,unsigned>::iterator,bool> insertProc;

  unsigned nEvts = 0;

  for (unsigned iF(0); iF<nFiles;++iF){

    std::ifstream input;
    std::ostringstream label;
    label << "/afs/cern.ch/work/a/amagnan/public/HGCalGeant4/git_V00-03-03/version_25/model_2/minbias/BON/run_0/momentum_list_layer30.dat";
    //label << "/afs/cern.ch/work/a/amagnan/public/HGCalGeant4/git_V00-03-03/version_25/model_2/minbias/BON/run_" << iF << "/layer30_version25_model2_BON_run" << iF << ".dat";
    input.open(label.str().c_str());
    if (!input.is_open()){
      std::cout << " Cannot open input file " << label.str() << ". Continue..." << std::endl;
      continue;
    }
    else std::cout << " Successfully opened file " << label.str() << std::endl;

    int pdgid=0;
    double time=0;
    double p=0;
    unsigned evtid=1000;
    unsigned trkS=0;
    double stepl=0;
    double totE=0;
    double nonIoE=0;
    std::string proc;
    double kinE=0;

    while (!input.eof()){
      std::string event("");
      std::getline(input,event);
      //std::cout << " -- line: " << event << std::endl;
      if (event.find("Event")!=event.npos){
	std::istringstream(event.substr(6,event.npos-6))>>evtid;
	//std::cout << " --Processing event id " << evtid 
	//<< " " << event.substr(6,event.npos-6)
	//<< std::endl;
	//return 1;
	nEvts++;
	continue;
      }
      std::istringstream(event)>>pdgid>>time>>p;
      if (pdgid>1000000000) std::istringstream(event)>>pdgid>>time>>p>>trkS>>stepl>>totE>>nonIoE>>proc>>kinE;
      //std::cout << evtid << " " << pdgid << " " << time << " " << p << std::endl;
      //return 1;
      if (evtid<1000 && pdgid!=0){
	if (p==0) {
	  if (abs(pdgid)>1000000000){
	    p_ionisingE->Fill(totE-nonIoE);
	    p_trkStatus->Fill(proc.c_str(),trkS);
	    p_stepLength->Fill(stepl);
	    insertProc = procnameMap.insert(std::pair<std::string,unsigned>(proc,1));
	    if (!insertProc.second) insertProc.first->second += 1;
	  }
	  continue;
	}
	isInserted = pdgIdMap.insert(std::pair<int,unsigned>(pdgid,1));
	if (!isInserted.second) {
	  //std::cout << " debug " << isInserted.first->first << " " << isInserted.first->second << std::endl;
	  isInserted.first->second += 1;
	}
	if (abs(pdgid)==2212) {
	  p_momentum_p->Fill(p/1000.);
	  p_logmomentum_p->Fill(log10(p/1000.));
	  p_time_p->Fill(time);
	  p_timevslogp_p->Fill(log10(p/1000.),time);
	  p_timevsp_p->Fill((p/1000.),time);
	}
	else if(abs(pdgid)==2112) {
	  p_momentum_n->Fill(p/1000.);
	  p_logmomentum_n->Fill(log10(p/1000.));
	  p_time_n->Fill(time);
	  p_timevslogp_n->Fill(log10(p/1000.),time);
	  p_timevsp_n->Fill((p/1000.),time);
	} 
	else if(abs(pdgid)==211) {
	  p_momentum_pi->Fill(p/1000.);
	  p_logmomentum_pi->Fill(log10(p/1000.));
	  p_time_pi->Fill(time);
	  p_timevslogp_pi->Fill(log10(p/1000.),time);
	  p_timevsp_pi->Fill((p/1000.),time);
	} 
	else if(abs(pdgid)>1000000000){// && abs(pdgid)<1000150000) {
	  p_momentum_Si->Fill(p/1000.);
	  p_logmomentum_Si->Fill(log10(p/1000.));
	  p_time_Si->Fill(time);
	  p_timevslogp_Si->Fill(log10(p/1000.),time);
	  p_timevsp_Si->Fill((p/1000.),time);
	  p_ionisingE->Fill(totE-nonIoE);
	  p_trkStatus->Fill(proc.c_str(),trkS);
	  p_stepLength->Fill(stepl);
	  insertProc = procnameMap.insert(std::pair<std::string,unsigned>(proc,1));
	  if (!insertProc.second) insertProc.first->second += 1;
	}
      }
      else {
	break;
      }
    }//read input file

  }//loop on files
  std::cout << " -- Read " << nEvts << " events." << std::endl;


  std::map<std::string,unsigned>::iterator iterProc = procnameMap.begin();
  for (;iterProc!=procnameMap.end();++iterProc){
    std::cout << iterProc->first << " " << iterProc->second << std::endl;
  }


  outputFile->cd();
  TGraph *gr = new TGraph();
  gr->SetName("grID");

  std::map<int,unsigned>::iterator iter = pdgIdMap.begin();
  unsigned i=0;
  unsigned nIons=0;
  std::vector<std::string> part;
  /*part.push_back("#bar{#Xi^{0}}");
  part.push_back("#Xi^{+}");
  part.push_back("#bar{#Sigma^{+}}");
  part.push_back("#bar{#Lambda}");
  part.push_back("#bar{#Sigma^{-}}");
  part.push_back("#bar{p}");
  part.push_back("#bar{n}");
  part.push_back("K^{-}");
  part.push_back("#pi^{-}");
  part.push_back("#pi^{0}");
  part.push_back("K^{0}_{L}");
  part.push_back("#pi^{+}");
  part.push_back("K^{0}_{S}");
  part.push_back("K^{+}");
  part.push_back("n");
  part.push_back("p");
  part.push_back("#Sigma^{-}");
  part.push_back("#Lambda");
  part.push_back("#Sigma^{+}");
  part.push_back("#Xi^{-}");
  part.push_back("#Xi^{0}");
  part.push_back("Ions");*/
  //part.push_back("#bar{#Lambda}");
  part.push_back("#bar{p}");
  part.push_back("#bar{n}");
  part.push_back("K^{-}");
  part.push_back("#pi^{-}");
  part.push_back("#bar{#nu_{#mu}}");
  part.push_back("#mu^{+}");
  part.push_back("#bar{#nu_{e}}");
  part.push_back("e^{+}");
  part.push_back("e^{-}");
  part.push_back("#nu_{e}");
  part.push_back("#mu^{-}");
  part.push_back("#nu_{#mu}");
  part.push_back("#gamma");
  part.push_back("K^{0}_{L}");
  part.push_back("#pi^{+}");
  part.push_back("K^{0}_{S}");
  part.push_back("K^{+}");
  part.push_back("n");
  part.push_back("p");
  part.push_back("#Sigma^{-}");
  part.push_back("#Lambda");
  part.push_back("Ions");
  for (;iter!=pdgIdMap.end();++iter){
    std::cout << i << " " << iter->first << " " << part[i] << " " << iter->second*1./nEvts << std::endl;
    //std::cout << i << " " << iter->first << " " << iter->second*1./nEvts << std::endl;
    if (iter->first<10000) {
      gr->SetPoint(i,i*1.+0.5,iter->second*1./nEvts);
      i++;
    }
    else {
      nIons += iter->second;
    }
  }
  gr->SetPoint(i,i*1.+0.5,nIons*1./nEvts);

  TAxis *ax = gr->GetHistogram()->GetXaxis();
  double x1 = ax->GetBinLowEdge(1);
  double x2 = ax->GetBinUpEdge(ax->GetNbins());
  gr->GetHistogram()->GetXaxis()->Set(gr->GetN(),x1,x2);
      
  for(unsigned k=0;k<gr->GetN();k++){
    gr->GetHistogram()->GetXaxis()->SetBinLabel(k+1,part[k].c_str());
  }

  outputFile->cd();
  gr->Write();
  outputFile->Write();
  return 0;

}//main
