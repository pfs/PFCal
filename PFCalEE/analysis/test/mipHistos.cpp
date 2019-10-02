#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<map>
#include <boost/algorithm/string.hpp>
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "boost/function.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TMath.h"
#include "TVectorD.h"

#include "HGCSSMipHit.hh"
#include "HGCSSGenParticle.hh"

#include "PositionFit.hh"
#include "SignalRegion.hh"

#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"

#include "utilities.h"

using boost::lexical_cast;
namespace po=boost::program_options;

int main(int argc, char** argv){//main  
  //Input output and config options
  std::string cfg;
  unsigned pNevts;
  std::string inFilePath;
  std::string outFilePath;
  unsigned nRunsIn;
  unsigned nRuns;
  unsigned startIdx;
  unsigned debug;
  bool do5lay;
  unsigned doTightSel;
  double etamean;
  double deta;

  po::options_description preconfig("Configuration"); 
  preconfig.add_options()("cfg,c",po::value<std::string>(&cfg)->required());
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(preconfig).allow_unregistered().run(), vm);
  po::notify(vm);
  po::options_description config("Configuration");
  config.add_options()
    //Input output and config options //->required()
    ("pNevts,n",       po::value<unsigned>(&pNevts)->default_value(0))
    ("inFilePath,i",   po::value<std::string>(&inFilePath)->required())
    ("outFilePath,o",  po::value<std::string>(&outFilePath)->required())
    ("nRunsIn",        po::value<unsigned>(&nRunsIn)->default_value(0))
    ("nRuns,R",        po::value<unsigned>(&nRuns)->default_value(0))
    ("startIdx,s",     po::value<unsigned>(&startIdx)->default_value(0))
    ("debug,d",        po::value<unsigned>(&debug)->default_value(0))
    ("do5lay",         po::value<bool>(&do5lay)->default_value(false))
    ("doTightSel",     po::value<unsigned>(&doTightSel)->default_value(0))
    ("etamean,e",      po::value<double>(&etamean)->default_value(2.0))
    ("deta",      po::value<double>(&deta)->default_value(0.1))
    ;

  po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
  po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
  po::notify(vm);


  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << inFilePath << std::endl
	    << " -- Output file path: " << outFilePath << std::endl
	    << " -- Do 5 layers: " << do5lay  << std::endl
	    << " -- Do tight selection (1: exactly one hit, 2: exactly in center of cells, 3: HCAL condition separation): " << doTightSel 
	    << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events in " << nRuns << "runs." << std::endl;
  else std::cout << pNevts << " events per run." << std::endl;

  const unsigned nNoise = 8;//5;//10;
  const double noise[nNoise] = {0,0.13,0.2,0.27,0.35,0.4,0.45,0.5};
  //for eta 2.7-3
  //const double noise[nNoise] = {0,0.13,0.2,0.27,0.35,0.4,0.5,0.6};
  //const double noise[nNoise] = {0,0.07,0.13,0.2,0.27,0.35,0.4,0.5,0.6,0.65};

  double threshBeforeAfterMin[nNoise];
  const double threshBeforeAfterMax = 2;
  std::cout << " ---- Selection settings: ---- " << std::endl
	    << " -------noise=threshBeforeAfterMin ";
  for (unsigned iN(0); iN<nNoise;++iN){
    threshBeforeAfterMin[iN] = std::max(0.9,2.2*noise[iN]);
    std::cout << noise[iN] << "=" << threshBeforeAfterMin[iN] << " ";
  }
  std::cout << std::endl
	    << " -------threshBeforeAfterMax " << threshBeforeAfterMax << std::endl
	    << " ------------------------------------------" << std::endl;



  const unsigned nLayers = 69;

  /////////////////////////////////////////////////////////////
  //input
  /////////////////////////////////////////////////////////////
  TChain *lTree = new TChain("MipTree");
  TFile * mipFile = 0;
  if (nRuns == 0){
    if (!testInputFile(inFilePath,mipFile)) return 1;
    lTree->AddFile(inFilePath.c_str());
  }
  else {
    for (unsigned i(startIdx);i<startIdx+nRuns;++i){
      std::ostringstream lstrrec;
      lstrrec << inFilePath << "_run" ;
      if (nRunsIn>0) lstrrec << nRunsIn*i << "_" << nRunsIn*(i+1)-1;
      else lstrrec << i ;
      lstrrec << ".root";
      if (!testInputFile(lstrrec.str(),mipFile)) continue;
      lTree->AddFile(lstrrec.str().c_str());
    }
  }
  if (!lTree){
    std::cout << " -- Error, tree MipTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  std::cout << " Trees added." << std::endl;


  TFile *outputFile = TFile::Open(outFilePath.c_str(),"RECREATE");
  
  if (!outputFile) {
    std::cout << " -- Error, output file " << outFilePath << " cannot be opened. Please create output directory. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
  }
  outputFile->cd();
  TTree *miptree=new TTree("MipTree","HGC MipAnalysis tree");

  HGCSSGenParticleVec lGenParts;
  miptree->Branch("HGCSSGenParticleVec","std::vector<HGCSSGenParticle>",&lGenParts);

  TH1F *hGenPnoSel = new TH1F("hGenPnoSel",";p^{gen} (GeV);Particles",5000,0,50);
  hGenPnoSel->StatOverflows(1);
  TH1F *hGenPTnoSel = new TH1F("hGenPTnoSel",";p_{T}^{gen} (GeV);Particles",2000,0,20);
  TH1F *hGenP = new TH1F("hGenP",";p^{gen} (GeV);Particles",5000,0,50000);
  hGenP->StatOverflows(1);
  TH1F *hGenPT = new TH1F("hGenPT",";p_{T}^{gen} (GeV);Particles",2000,0,20000);
  TH1F *hGenEta = new TH1F("hGenEta",";#eta^{gen};Particles",32,1.4,3);
  TH1F *hGenPhi = new TH1F("hGenPhi",";#phi^{gen};Particles",180,-3.1416,3.1416);
  TH1F *hGenPID = new TH1F("hGenPID",";pdgId^{gen};Particles",1000,-500,500);

  TH1F *hE[nNoise][nLayers+1];
  TH1F *hEsig[nNoise][nLayers+1];
  TH1F *hEmaxNeigh[nNoise][nLayers+1];
  TH1F *hEmaxNeighsig[nNoise][nLayers+1];
  TH1F *hEsumNeigh[nNoise][nLayers+1];
  TH1F *hEsumNeighsig[nNoise][nLayers+1];
  TH1F *hnTrkPerEvt[nNoise][nLayers+1];
  TH1F *hNeighIdx[nNoise][nLayers+1][4];
  TH1F *hNeighIdxsig[nNoise][nLayers+1][4];
  TH1F *hNeighE[nNoise][nLayers+1][4];
  TH1F *hNeighEsig[nNoise][nLayers+1][4];

  for (unsigned iN(0); iN<nNoise;++iN){
    for (unsigned iL(0); iL<nLayers+1; ++iL){
      std::ostringstream label;
      label << "Ehit_" << iN << "_" << iL;
      hE[iN][iL] = new TH1F(label.str().c_str(),";E (mips); Mip tracks",100,0.5,5);
      label.str("");
      label << "Ehitsig_" << iN << "_" << iL;
      hEsig[iN][iL] = new TH1F(label.str().c_str(),";E (mips); Mip tracks",100,0.5,5);
      label.str("");
      label << "EmaxNeigh_" << iN << "_" << iL;
      hEmaxNeigh[iN][iL] = new TH1F(label.str().c_str(),";E_{max}^{neigh} (mips); Mip tracks",100,0,5);
      label.str("");
      label << "EmaxNeighsig_" << iN << "_" << iL;
      hEmaxNeighsig[iN][iL] = new TH1F(label.str().c_str(),";E_{max}^{neigh} (mips); Mip tracks",100,0,5);
      label.str("");
      label << "EsumNeigh_" << iN << "_" << iL;
      hEsumNeigh[iN][iL] = new TH1F(label.str().c_str(),";E_{sum}^{neigh} (mips); Mip tracks",100,0,5);
      label.str("");
      label << "EsumNeighsig_" << iN << "_" << iL;
      hEsumNeighsig[iN][iL] = new TH1F(label.str().c_str(),";E_{sum}^{neigh} (mips); Mip tracks",100,0,5);
      label.str("");
      label << "nTrkPerEvt_" << iN << "_" << iL;
      hnTrkPerEvt[iN][iL] = new TH1F(label.str().c_str(),";n_{trks}/event; Events",iL==nLayers?2000:100,0,iL==nLayers?2000:100);
      //index of neighbour in range for L-2->L+2
      for (unsigned jL(0);jL<4;++jL){
	label.str("");
	label << "neighIdx_" << iN << "_" << iL << "_neighLay" << jL;
	hNeighIdx[iN][iL][jL] = new TH1F(label.str().c_str(),";Neighbour Idx; Mip tracks",9,0,9);
	label.str("");
	label << "neighE_" << iN << "_" << iL << "_neighLay" << jL;
	hNeighE[iN][iL][jL] = new TH1F(label.str().c_str(),";Neighbour E (mips); Mip tracks",100,0,5);
	label.str("");
	label << "neighIdxsig_" << iN << "_" << iL << "_neighLay" << jL;
	hNeighIdxsig[iN][iL][jL] = new TH1F(label.str().c_str(),";Neighbour Idx; Mip tracks",9,0,9);
	label.str("");
	label << "neighEsig_" << iN << "_" << iL << "_neighLay" << jL;
	hNeighEsig[iN][iL][jL] = new TH1F(label.str().c_str(),";Neighbour E (mips); Mip tracks",100,0,5);
       }
    }
  }

  std::vector<HGCSSMipHitVec*> miphitvec;
  miphitvec.resize(nNoise,0);
  std::vector<HGCSSGenParticle> * genvec = 0;

  if (lTree->GetBranch("HGCSSGenParticleVec")) lTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);


  for (unsigned iN(0); iN<nNoise;++iN){
    std::ostringstream label;

    label << "MipHitVec_noise" << iN;//noise[iN];
    lTree->SetBranchAddress(label.str().c_str(),&miphitvec[iN]);
  }
  
  const unsigned nEvts = ((pNevts > lTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lTree->GetEntries()) : pNevts) ;
  
  std::cout << " -- Processing " << nEvts << " events out of " << lTree->GetEntries() << std::endl;
  bool firstprint = true;
  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
    lTree->GetEntry(ievt);
    //save genparticles
    lGenParts.clear();
    unsigned nGenParts = (*genvec).size();
    lGenParts.reserve(nGenParts);
    //if (nGenParts != 1) continue;
    for (unsigned iP(0); iP<nGenParts; ++iP){//loop on gen particles
      HGCSSGenParticle lPart = (*genvec)[iP];
      lGenParts.push_back(lPart);
      hGenPnoSel->Fill(lPart.p()/1000.);
      hGenPTnoSel->Fill(lPart.pt()/1000.);
    }

    miptree->Fill();

    bool passCut = false;
    for (unsigned iN(0); iN<nNoise;++iN){//loop on noise
      unsigned nTrkPerEvt[nLayers+1];
      for (unsigned iL(0); iL<nLayers+1; ++iL){
	nTrkPerEvt[iL] = 0;
      }
      if (debug) std::cout << " -- noise " << iN << " nGenParts = " << nGenParts << " nHits = " << (*miphitvec[iN]).size() << std::endl;
      for (unsigned iH(0); iH<(*miphitvec[iN]).size(); ++iH){//loop on hits
	HGCSSMipHit lHit = (*miphitvec[iN])[iH];
	unsigned layer = lHit.l();
	//fill histos
	bool fail = false;
	//count neighbours within range
	unsigned nNeigh[4] = {0,0,0,0};
	unsigned neiId[4] = {0,0,0,0};
	double neiE[4] = {0,0,0,0};
	unsigned nNei = layer<53?7:9;
	if (doTightSel<3 || layer>26){
	  for (unsigned ineigh(0); ineigh<nNei; ++ineigh){
	    if (lHit.neigh_e_prev2layer(ineigh)>threshBeforeAfterMin[iN] && lHit.neigh_e_prev2layer(ineigh)<threshBeforeAfterMax) {
	      nNeigh[0]++;
	      neiId[0] = ineigh;
	      neiE[0] = lHit.neigh_e_prev2layer(ineigh);
	    }
	    if (lHit.neigh_e_prevlayer(ineigh)>threshBeforeAfterMin[iN] && lHit.neigh_e_prevlayer(ineigh)<threshBeforeAfterMax) {
	      nNeigh[1]++;
	      neiId[1] = ineigh;
	      neiE[1] = lHit.neigh_e_prevlayer(ineigh);
	    }
	    if (lHit.neigh_e_nextlayer(ineigh)>threshBeforeAfterMin[iN] && lHit.neigh_e_nextlayer(ineigh)<threshBeforeAfterMax) {
	      nNeigh[2]++;
	      neiId[2] = ineigh;
	      neiE[2] = lHit.neigh_e_nextlayer(ineigh);
	    }
	    if (lHit.neigh_e_next2layer(ineigh)>threshBeforeAfterMin[iN] && lHit.neigh_e_next2layer(ineigh)<threshBeforeAfterMax) {
	      nNeigh[3]++;
	      neiId[3] = ineigh;
	      neiE[3] = lHit.neigh_e_next2layer(ineigh);
	    }
	  }
	}
	else {//hcal neigh distance selection
	  for (unsigned jH(0); jH<(*miphitvec[iN]).size(); ++jH){//loop on hits
	    if (jH==iH) continue;
	    HGCSSMipHit lSecHit = (*miphitvec[iN])[jH];
	    unsigned seclay = lSecHit.l();
	    if (seclay == layer || layer>26) continue;
	    bool passLay = 
	      (layer>2 && layer<25 && abs(layer-seclay)==3) ||
	      (layer==2 && ((layer-seclay)==2 || (seclay-layer)==3)) ||
	      (layer<2 && (seclay-layer)==3) ||
	      (layer==25 && ((seclay-layer)==2 || (layer-seclay)==3))||
	      (layer>25 && layer-seclay==3);
	    if (!passLay) continue;
	    double secE = lSecHit.e();
	    ROOT::Math::XYZPoint pos = ROOT::Math::XYZPoint(lHit.x(),lHit.y(),lHit.z());
	    ROOT::Math::XYZPoint secpos = ROOT::Math::XYZPoint(lSecHit.x(),lSecHit.y(),lSecHit.z());
	    if (secE>threshBeforeAfterMin[iN] && secE<threshBeforeAfterMax && fabs(pos.eta()-secpos.eta()) < 0.01){
	      if ( int(layer-seclay) < 0 ){
		nNeigh[2]++;
		neiE[2] = secE;
		neiId[2] = 0;
	      }
	      else {
		nNeigh[1]++;
		neiId[1] = 0;
		neiE[1] = secE;
	      }
	    }//energy in correct range
	  }//loop on hits
	}//HCAL hit selection
	if (do5lay){
	  if (!doTightSel) fail = lHit.getMaxEnergy(-2)<threshBeforeAfterMin[iN] || lHit.getMaxEnergy(2)<threshBeforeAfterMin[iN] || lHit.getMaxEnergy(-2)>threshBeforeAfterMax || lHit.getMaxEnergy(2)>threshBeforeAfterMax;
	  else {//doTightSel

	    if (doTightSel==1) {
	      bool fail0 = nNeigh[0]!=1 || nNeigh[1]!=1 ||nNeigh[2]!=1 || nNeigh[3]!=1;
	      if (layer<53){
		bool pass0 = neiId[0]==0 && neiId[1]==0 && neiId[2]==0 && neiId[3]==0 ;//((neiId[2]==0 && (neiId[3]==0||neiId[3]>=4))  || (neiId[2]>=4 && neiId[3]==neiId[2]));
		bool pass1 = neiId[0]==1 && (neiId[1]==1 || neiId[1]==0) && ((neiId[2]==0 && (neiId[3]==0||neiId[3]==6)) || (neiId[2]==6 && neiId[3]==6)) ;
		bool pass2 = neiId[0]==2 && (neiId[1]==2 || neiId[1]==0) && ((neiId[2]==0 && (neiId[3]==0||neiId[3]==5)) || (neiId[2]==5 && neiId[3]==5)) ;
		bool pass3 = neiId[0]==3 && (neiId[1]==3 || neiId[1]==0) && ((neiId[2]==0 && (neiId[3]==0||neiId[3]==4)) || (neiId[2]==4 && neiId[3]==4)) ;
		bool pass4 = neiId[0]==4 && (neiId[1]==4 || neiId[1]==0) && ((neiId[2]==0 && (neiId[3]==0||neiId[3]==3)) || (neiId[2]==3 && neiId[3]==3)) ;
		bool pass5 = neiId[0]==5 && (neiId[1]==5 || neiId[1]==0) && ((neiId[2]==0 && (neiId[3]==0||neiId[3]==2)) || (neiId[2]==2 && neiId[3]==2)) ;
		bool pass6 = neiId[0]==6 && (neiId[1]==6 || neiId[1]==0) && ((neiId[2]==0 && (neiId[3]==0||neiId[3]==1)) || (neiId[2]==1 && neiId[3]==1)) ;
		fail = fail0 || !(pass0||pass1||pass2||pass3||pass4||pass5||pass6);
	      //fail = nNeigh[0]!=1 || nNeigh[1]!=1 ||nNeigh[2]!=1 ||nNeigh[3]!=1;
	      }
	      else {
		bool pass0 = neiId[0]==0 && neiId[1]==0 && neiId[2]==0 && neiId[3]==0 ;//((neiId[2]==0 && (neiId[3]==0||neiId[3]>=4))  || (neiId[2]>=4 && neiId[3]==neiId[2]));
		bool pass1 = neiId[0]==1 && (neiId[1]==1 || neiId[1]==0) && ((neiId[2]==0 && (neiId[3]==0||neiId[3]==8)) || (neiId[2]==8 && neiId[3]==8)) ;
		bool pass2 = neiId[0]==2 && (neiId[1]==2 || neiId[1]==0) && ((neiId[2]==0 && (neiId[3]==0||neiId[3]==7)) || (neiId[2]==7 && neiId[3]==7)) ;
		bool pass3 = neiId[0]==3 && (neiId[1]==3 || neiId[1]==0) && ((neiId[2]==0 && (neiId[3]==0||neiId[3]==6)) || (neiId[2]==6 && neiId[3]==6)) ;
		bool pass4 = neiId[0]==4 && (neiId[1]==4 || neiId[1]==0) && ((neiId[2]==0 && (neiId[3]==0||neiId[3]==5)) || (neiId[2]==5 && neiId[3]==5)) ;
		bool pass5 = neiId[0]==5 && (neiId[1]==5 || neiId[1]==0) && ((neiId[2]==0 && (neiId[3]==0||neiId[3]==4)) || (neiId[2]==4 && neiId[3]==4)) ;
		bool pass6 = neiId[0]==6 && (neiId[1]==6 || neiId[1]==0) && ((neiId[2]==0 && (neiId[3]==0||neiId[3]==3)) || (neiId[2]==3 && neiId[3]==3)) ;
		bool pass7 = neiId[0]==7 && (neiId[1]==7 || neiId[1]==0) && ((neiId[2]==0 && (neiId[3]==0||neiId[3]==2)) || (neiId[2]==2 && neiId[3]==2)) ;
		bool pass8 = neiId[0]==8 && (neiId[1]==8 || neiId[1]==0) && ((neiId[2]==0 && (neiId[3]==0||neiId[3]==1)) || (neiId[2]==1 && neiId[3]==1)) ;
		fail = fail0 || !(pass0||pass1||pass2||pass3||pass4||pass5||pass6||pass7||pass8);

	      }
	      
	    }//tightsel
	    else if (doTightSel==2) {//tightsel==2: want exactly the layer aligned and nothing else...
	      bool failneigh = false;
	      //fail if one of the neighbour in all 4 layers around is within range
	      for (unsigned ineigh(0); ineigh<nNei-1; ++ineigh){
		bool tmp = 
		  (lHit.neigh_e_prev2layer(ineigh+1)>threshBeforeAfterMin[iN] && lHit.neigh_e_prev2layer(ineigh+1)<threshBeforeAfterMax) || 
		  (lHit.neigh_e_prevlayer(ineigh+1)>threshBeforeAfterMin[iN] && lHit.neigh_e_prevlayer(ineigh+1)<threshBeforeAfterMax) || 
		  (lHit.neigh_e_next2layer(ineigh+1)>threshBeforeAfterMin[iN] && lHit.neigh_e_next2layer(ineigh+1)<threshBeforeAfterMax) || 
		  (lHit.neigh_e_nextlayer(ineigh+1)>threshBeforeAfterMin[iN] && lHit.neigh_e_nextlayer(ineigh+1)<threshBeforeAfterMax);
		if (tmp) {
		  failneigh = true;
		  break;
		}
	      }
	      //fail if hit aligned with layer in all 4 layers around is not in range
	      bool centralfail = lHit.neigh_e_prev2layer(0)<threshBeforeAfterMin[iN] || lHit.neigh_e_next2layer(0)<threshBeforeAfterMin[iN] || lHit.neigh_e_prev2layer(0)>threshBeforeAfterMax || lHit.neigh_e_next2layer(0)>threshBeforeAfterMax || lHit.neigh_e_prevlayer(0)<threshBeforeAfterMin[iN] || lHit.neigh_e_nextlayer(0)<threshBeforeAfterMin[iN] || lHit.neigh_e_prevlayer(0)>threshBeforeAfterMax || lHit.neigh_e_nextlayer(0)>threshBeforeAfterMax ;
	      fail = failneigh || centralfail;
	    }
	  }//doTightSel
	}//do5lay
	else {
	  if (doTightSel){
	    if (doTightSel==1 || doTightSel==3) {
	      if (debug>1) std::cout << " -- Hit #" << iH << " Number of neighbours " <<  nNeigh[1] << " " <<  nNeigh[2] << std::endl;
	      bool fail0 = nNeigh[1]!=1 || nNeigh[2]!=1;
	      if (layer<53){
		bool pass0 = neiId[1]==0 && neiId[2]==0;
		bool pass1 = neiId[1]==1 && (neiId[2]==0 || neiId[2]==6);
		bool pass2 = neiId[1]==2 && (neiId[2]==0 || neiId[2]==5);
		bool pass3 = neiId[1]==3 && (neiId[2]==0 || neiId[2]==4);
		bool pass4 = neiId[1]==4 && (neiId[2]==0 || neiId[2]==3);
		bool pass5 = neiId[1]==5 && (neiId[2]==0 || neiId[2]==2);
		bool pass6 = neiId[1]==6 && (neiId[2]==0 || neiId[2]==1);
		fail = fail0 || !(pass0||pass1||pass2||pass3||pass4||pass5||pass6);
	      }
	      else {
		bool pass0 = neiId[1]==0 && neiId[2]==0;
		bool pass1 = neiId[1]==1 && (neiId[2]==0 || neiId[2]==8);
		bool pass2 = neiId[1]==2 && (neiId[2]==0 || neiId[2]==7);
		bool pass3 = neiId[1]==3 && (neiId[2]==0 || neiId[2]==6);
		bool pass4 = neiId[1]==4 && (neiId[2]==0 || neiId[2]==5);
		bool pass5 = neiId[1]==5 && (neiId[2]==0 || neiId[2]==4);
		bool pass6 = neiId[1]==6 && (neiId[2]==0 || neiId[2]==3);
		bool pass7 = neiId[1]==7 && (neiId[2]==0 || neiId[2]==2);
		bool pass8 = neiId[1]==8 && (neiId[2]==0 || neiId[2]==1);
		fail = fail0 || !(pass0||pass1||pass2||pass3||pass4||pass5||pass6||pass7||pass8);
	      }
	    }
	    else if (doTightSel==2){
	      bool failneigh = false;
	      //fail if one of the neighbour in all 2 layers around is within range
	      for (unsigned ineigh(0); ineigh<nNei-1; ++ineigh){
		bool tmp = 
		  (lHit.neigh_e_prevlayer(ineigh+1)>threshBeforeAfterMin[iN] && lHit.neigh_e_prevlayer(ineigh+1)<threshBeforeAfterMax) || 
		  (lHit.neigh_e_nextlayer(ineigh+1)>threshBeforeAfterMin[iN] && lHit.neigh_e_nextlayer(ineigh+1)<threshBeforeAfterMax);
		if (tmp) {
		  failneigh = true;
		  break;
		}
	      }
	      bool centralfail = lHit.neigh_e_prevlayer(0)<threshBeforeAfterMin[iN] || lHit.neigh_e_nextlayer(0)<threshBeforeAfterMin[iN] || lHit.neigh_e_prevlayer(0)>threshBeforeAfterMax || lHit.neigh_e_nextlayer(0)>threshBeforeAfterMax ;
	      fail = failneigh || centralfail;
	    }
	  }//doTightSel
	}//3 lay
	if (fail) continue;

	if (iN==0) passCut = true;

	ROOT::Math::XYZPoint pos = ROOT::Math::XYZPoint(lHit.x(),lHit.y(),lHit.z());
	if (iN==0 && debug) std::cout << " - layer " << layer << " sel hit eta = " << pos.eta() << " phi = " << pos.phi() << std::endl;

	hE[iN][layer]->Fill(lHit.e());
	hE[iN][nLayers]->Fill(lHit.e());
	hEmaxNeigh[iN][layer]->Fill(lHit.getMaxEnergy(0));
	hEmaxNeigh[iN][nLayers]->Fill(lHit.getMaxEnergy(0));
	hEsumNeigh[iN][layer]->Fill(lHit.getSumEnergy(0));
	hEsumNeigh[iN][nLayers]->Fill(lHit.getSumEnergy(0));

	for (unsigned jL(0);jL<4;++jL){
	  if (nNeigh[jL]!=1) continue;
	  hNeighIdx[iN][layer][jL]->Fill(neiId[jL]);
	  hNeighE[iN][layer][jL]->Fill(neiE[jL]);
	  hNeighIdx[iN][nLayers][jL]->Fill(neiId[jL]);
	  hNeighE[iN][nLayers][jL]->Fill(neiE[jL]);
	}
	if (lHit.noiseFrac()<0.5) {
	  hEsig[iN][layer]->Fill(lHit.e());
	  hEsig[iN][nLayers]->Fill(lHit.e());
	  hEmaxNeighsig[iN][layer]->Fill(lHit.getMaxEnergy(0));
	  hEmaxNeighsig[iN][nLayers]->Fill(lHit.getMaxEnergy(0));
	  hEsumNeighsig[iN][layer]->Fill(lHit.getSumEnergy(0));
	  hEsumNeighsig[iN][nLayers]->Fill(lHit.getSumEnergy(0));
	  for (unsigned jL(0);jL<4;++jL){
	    if (nNeigh[jL]!=1) continue;
	    hNeighIdxsig[iN][layer][jL]->Fill(neiId[jL]);
	    hNeighEsig[iN][layer][jL]->Fill(neiE[jL]);
	    hNeighIdxsig[iN][nLayers][jL]->Fill(neiId[jL]);
	    hNeighEsig[iN][nLayers][jL]->Fill(neiE[jL]);
	  }
	} 
	/* else {
	   if (firstprint && layer==6) {
	   std::cout << " -- Event " << ievt << " noise " << noise[iN] << std::endl
		      << "double neighE_prev2[7] = {";
	    for (unsigned ineigh(0); ineigh<7; ++ineigh){
	      std::cout << lHit.neigh_e_prev2layer(ineigh);
	      if (ineigh<6) std::cout << ",";
	    }
	    std::cout << "};" << std::endl
		      << "double neighE_prev[7] = {";
	    for (unsigned ineigh(0); ineigh<7; ++ineigh){
	      std::cout << lHit.neigh_e_prevlayer(ineigh);
	      if (ineigh<6) std::cout << ",";
	    }
	    std::cout << "};" << std::endl
		      << "double neighE_same[7] = {" << lHit.e() << ",";
	    for (unsigned ineigh(1); ineigh<7; ++ineigh){
	      std::cout << lHit.neigh_e_samelayer(ineigh);
	      if (ineigh<6) std::cout << ",";
	    }
	    std::cout << "};" << std::endl
		      << "double neighE_next[7] = {";
	    for (unsigned ineigh(0); ineigh<7; ++ineigh){
	      std::cout << lHit.neigh_e_nextlayer(ineigh);
	      if (ineigh<6) std::cout << ",";
	    }
	    std::cout << "};" << std::endl
		      << "double neighE_next2[7] = {";
	    for (unsigned ineigh(0); ineigh<7; ++ineigh){
	      std::cout << lHit.neigh_e_next2layer(ineigh);
	      if (ineigh<6) std::cout << ",";
	    }
	    std::cout << "};" << std::endl;
	    if (iN==nNoise-1) firstprint = false;
	  }
	  
	}*/

	nTrkPerEvt[layer]++;
	nTrkPerEvt[nLayers]++;
      }//loop on hits


      for (unsigned iL(0); iL<nLayers+1; ++iL){
	hnTrkPerEvt[iN][iL]->Fill(nTrkPerEvt[iL]);
      }
    }//loop on noise

    //find genparticles corresponding to selected hits...
    if (passCut){
      for (unsigned iP(0); iP<(*genvec).size(); ++iP){//loop on gen particles
	HGCSSGenParticle lPart = (*genvec)[iP];
	hGenP->Fill(lPart.p()/1000.);
	hGenPT->Fill(lPart.pt()/1000.);
	hGenEta->Fill(lPart.eta());
	hGenPhi->Fill(lPart.phi());
	hGenPID->Fill(lPart.pdgid());
	if (debug) std::cout << " --- sel genParticle eta = " << lPart.eta() << " phi = " << lPart.phi() << " E = " << lPart.E()/1000. << " GeV" << std::endl;
	
      }
    }

  }//loop on entries

  outputFile->cd();
  miptree->Write();
  outputFile->Write();
  outputFile->Close();
  
  return 0;


}//main
