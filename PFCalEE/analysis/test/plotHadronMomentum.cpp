#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TStyle.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"
#include <boost/algorithm/string.hpp>
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "boost/function.hpp"

#include "HGCSSEvent.hh"
#include "HGCSSInfo.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSParameters.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDigitisation.hh"
#include "HGCSSDetector.hh"
#include "HGCSSGeometryConversion.hh"
#include "HGCSSPUenergy.hh"

#include "utilities.h"

using boost::lexical_cast;
namespace po=boost::program_options;


int main(int argc, char** argv){//main
  //Input output and config options
  std::string cfg;
  unsigned pNevts;
  std::string filePath;
  std::string digifilePath;
  unsigned nRuns;
  std::string simFileName;
  std::string recoFileName;
  std::string outPath;
  unsigned debug;
  unsigned start;
  po::options_description preconfig("Configuration"); 
  preconfig.add_options()("cfg,c",po::value<std::string>(&cfg)->required());
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(preconfig).allow_unregistered().run(), vm);
  po::notify(vm);
  po::options_description config("Configuration");
  config.add_options()
    ("pNevts,n",       po::value<unsigned>(&pNevts)->default_value(0))
    ("filePath,i",     po::value<std::string>(&filePath)->required())
    ("digifilePath", po::value<std::string>(&digifilePath)->default_value(""))
    ("nRuns",        po::value<unsigned>(&nRuns)->default_value(0))
    ("simFileName,s",  po::value<std::string>(&simFileName)->required())
    ("recoFileName,r", po::value<std::string>(&recoFileName)->required())
    ("outPath,o",      po::value<std::string>(&outPath)->required())
    ("debug,d",        po::value<unsigned>(&debug)->default_value(0))
    ("start",          po::value<unsigned>(&start)->default_value(0))
    ;

  po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
  po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
  po::notify(vm);



  std::string inFilePath = filePath+simFileName;
  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << filePath << std::endl
	    << " -- Digi Input file path: " << digifilePath << std::endl
	    << " -- Output file path: " << outPath << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events." << std::endl;
  /////////////////////////////////////////////////////////////
  //input
  /////////////////////////////////////////////////////////////

  std::ostringstream inputsim;
  inputsim << filePath << "/" << simFileName;
  std::ostringstream inputrec;
  if (digifilePath.size()==0)
    inputrec << filePath << "/" << recoFileName;
  else 
    inputrec << digifilePath << "/" << recoFileName;

  //std::cout << inputsim.str() << " " << inputrec.str() << std::endl;

  HGCSSInfo * info;

  TChain *lSimTree = new TChain("HGCSSTree");
  TChain *lRecTree = 0;
  
  TFile * simFile = 0;
  TFile * recFile = 0;

  if (recoFileName.find("Digi") != recoFileName.npos) 
    lRecTree = new TChain("RecoTree");
  else lRecTree = new TChain("PUTree");

  if (nRuns == 0){
    if (!testInputFile(inputsim.str(),simFile)) return 1;
    lSimTree->AddFile(inputsim.str().c_str());
    if (!testInputFile(inputrec.str(),recFile)) return 1;
    lRecTree->AddFile(inputrec.str().c_str());
  }
  else {
    for (unsigned i(start);i<start+nRuns;++i){
      std::ostringstream lstr;
      lstr << inputsim.str() << "_" << i << ".root";
      if (testInputFile(lstr.str(),simFile)){  
	if (simFile) info =(HGCSSInfo*)simFile->Get("Info");
	else {
	  std::cout << " -- Error in getting information from simfile!" << std::endl;
	  return 1;
	}
      }
      else continue;
      lSimTree->AddFile(lstr.str().c_str());
      lstr.str("");
      lstr << inputrec.str() << "_" << i << ".root";
      if (!testInputFile(lstr.str(),recFile)) continue;
      lRecTree->AddFile(lstr.str().c_str());
    }
  }

  if (!lSimTree){
    std::cout << " -- Error, tree HGCSSTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  if (!lRecTree){
    std::cout << " -- Error, tree RecoTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }


  /////////////////////////////////////////////////////////////
  //Info
  /////////////////////////////////////////////////////////////

  const double cellSize = info->cellSize();
  const unsigned versionNumber = info->version();
  const unsigned model = info->model();
  
  std::cout << " -- Version number is : " << versionNumber 
	    << ", model = " << model
	    << ", cellSize = " << cellSize
	    << std::endl;


  //initialise detector
  HGCSSDetector & myDetector = theDetector();
 
  myDetector.buildDetector(versionNumber,true,false);

  const unsigned nLayers = myDetector.nLayers();
  const unsigned nSections = myDetector.nSections();

  std::cout << " -- N layers = " << nLayers << std::endl
	    << " -- N sections = " << nSections << std::endl;

  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  ///////// Output File // /////////////////////////
  //////////////////////////////////////////////////
  //////////////////////////////////////////////////


  std::ostringstream lstr;
  lstr << outPath << "_" << start << ".root";
  TFile *outputFile = TFile::Open(lstr.str().c_str(),"RECREATE");
  
  if (!outputFile) {
    std::cout << " -- Error, output file " << outPath << " cannot be opened. Please create output directory. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
  }
  outputFile->cd();

  outputFile->mkdir("pipm");
  outputFile->mkdir("protons");
  outputFile->mkdir("neutrons");
  outputFile->cd("protons");
  TH1F *p_gFrac_p = new TH1F("p_gFrac_p",";#gamma fraction;hits",101,0,1.01);
  TH1F *p_eFrac_p = new TH1F("p_eFrac_p",";e fraction;hits",101,0,1.01);
  TH1F *p_muFrac_p = new TH1F("p_muFrac_p",";#mu fraction;hits",101,0,1.01);
  TH1F *p_hadFrac_p = new TH1F("p_hadFrac_p",";had fraction;hits",101,0,1.01);
  TH1F *p_neutronFrac_p = new TH1F("p_neutronFrac_p",";neutron fraction;hits",101,0,1.01);
  TH1F *p_protonFrac_p = new TH1F("p_protonFrac_p",";proton fraction;hits",101,0,1.01);
  TH1F *p_logmomentum_p  = new TH1F("p_logmomentum_p",";log(p) (log(GeV));particles",100,-2,2);
  TH1F *p_momentum_p  = new TH1F("p_momentum_p",";p (GeV);particles",1000,0,100);
  TH1F *p_mainParentEfrac_p  = new TH1F("p_mainParentEfrac_p",";eParent/etot;particles",100,0,10);

  outputFile->cd();
  outputFile->cd("neutrons");
  TH1F *p_gFrac_n = new TH1F("p_gFrac_n",";#gamma fraction;hits",101,0,1.01);
  TH1F *p_eFrac_n = new TH1F("p_eFrac_n",";e fraction;hits",101,0,1.01);
  TH1F *p_muFrac_n = new TH1F("p_muFrac_n",";#mu fraction;hits",101,0,1.01);
  TH1F *p_hadFrac_n = new TH1F("p_hadFrac_n",";had fraction;hits",101,0,1.01);
  TH1F *p_neutronFrac_n = new TH1F("p_neutronFrac_n",";neutron fraction;hits",101,0,1.01);
  TH1F *p_protonFrac_n = new TH1F("p_protonFrac_n",";proton fraction;hits",101,0,1.01);
  TH1F *p_logmomentum_n  = new TH1F("p_logmomentum_n",";log(p) (log(GeV));particles",100,-2,2);
  TH1F *p_momentum_n  = new TH1F("p_momentum_n",";p (GeV);particles",1000,0,100);
  TH1F *p_mainParentEfrac_n  = new TH1F("p_mainParentEfrac_n",";eParent/etot;particles",100,0,10);

  outputFile->cd();
  outputFile->cd("pipm");
  TH1F *p_gFrac_pi = new TH1F("p_gFrac_pi",";#gamma fraction;hits",101,0,1.01);
  TH1F *p_eFrac_pi = new TH1F("p_eFrac_pi",";e fraction;hits",101,0,1.01);
  TH1F *p_muFrac_pi = new TH1F("p_muFrac_pi",";#mu fraction;hits",101,0,1.01);
  TH1F *p_hadFrac_pi = new TH1F("p_hadFrac_pi",";had fraction;hits",101,0,1.01);
  TH1F *p_neutronFrac_pi = new TH1F("p_neutronFrac_pi",";neutron fraction;hits",101,0,1.01);
  TH1F *p_protonFrac_pi = new TH1F("p_protonFrac_pi",";proton fraction;hits",101,0,1.01);
  TH1F *p_logmomentum_pi  = new TH1F("p_logmomentum_pi",";log(p) (log(GeV));particles",100,-2,2);
  TH1F *p_momentum_pi  = new TH1F("p_momentum_pi",";p (GeV);particles",1000,0,100);
  TH1F *p_mainParentEfrac_pi  = new TH1F("p_mainParentEfrac_pi",";eParent/etot;particles",100,0,10);

 
  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  //////////////////////////////////////////////////


  const unsigned nEvts = ((pNevts > lSimTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lSimTree->GetEntries()) : pNevts) ;
    //loop on events
  HGCSSEvent * event = 0;
  std::vector<HGCSSSamplingSection> * ssvec = 0;
  std::vector<HGCSSSimHit> * simhitvec = 0;
  std::vector<HGCSSRecoHit> * rechitvec = 0;
  std::vector<HGCSSGenParticle> * genvec = 0;
  unsigned nPuVtx = 0;

  lSimTree->SetBranchAddress("HGCSSEvent",&event);
  lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
  lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
  lSimTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);

  lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  if (lRecTree->GetBranch("nPuVtx")) lRecTree->SetBranchAddress("nPuVtx",&nPuVtx);

  unsigned nProtons = 0;
  unsigned nPipm = 0;
  unsigned nNeutrons = 0;

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;

    lSimTree->GetEntry(ievt);
    for (unsigned iG(0); iG<(*genvec).size(); ++iG){//loop on genparticles
      const HGCSSGenParticle & lgen = (*genvec)[iG];
      if (lgen.eta()<2.8 || lgen.eta()>3.0) continue;
      unsigned pdgid = abs(lgen.pdgid());
      if (pdgid != 2212 && pdgid != 2112 && pdgid != 211) continue;
      unsigned gentrkId =  lgen.trackID();
      for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
	const HGCSSSimHit & lHit = (*simhitvec)[iH];
	if (lHit.layer()!=30) continue;
	unsigned trkId =  lHit.mainParentTrackID();
	if (gentrkId != trkId) continue;
	double p = sqrt(pow(lgen.px(),2)+pow(lgen.py(),2)+pow(lgen.pz(),2));
	if (pdgid == 211) {
	  p_logmomentum_pi->Fill(log10(p/1000.));
	  p_momentum_pi->Fill(p/1000.);
	  if (lHit.energy()>0) p_mainParentEfrac_pi->Fill(lHit.mainParentEfrac());
	  p_protonFrac_pi->Fill(lHit.protonFrac());
	  p_neutronFrac_pi->Fill(lHit.neutronFrac());
	  p_gFrac_pi->Fill(lHit.gFrac());
	  p_eFrac_pi->Fill(lHit.eFrac());
	  p_muFrac_pi->Fill(lHit.muFrac());
	  p_hadFrac_pi->Fill(lHit.hadFrac());
	  if (debug){
	    std::cout << " -- Found proton! ";
	    lHit.Print(std::cout);
	  }
	  nPipm++;
	  break;
	}
	else if (pdgid == 2212) {
	  p_logmomentum_p->Fill(log10(p/1000.));
	  p_momentum_p->Fill(p/1000.);
	  if (lHit.energy()>0) p_mainParentEfrac_p->Fill(lHit.mainParentEfrac());
	  p_protonFrac_p->Fill(lHit.protonFrac());
	  p_neutronFrac_p->Fill(lHit.neutronFrac());
	  p_gFrac_p->Fill(lHit.gFrac());
	  p_eFrac_p->Fill(lHit.eFrac());
	  p_muFrac_p->Fill(lHit.muFrac());
	  p_hadFrac_p->Fill(lHit.hadFrac());
	  if (debug){
	    std::cout << " -- Found proton! ";
	    lHit.Print(std::cout);
	  }
	  nProtons++;
	  break;
	}
	else if (pdgid == 2112) {
	  p_logmomentum_n->Fill(log10(p/1000.));
	  p_momentum_n->Fill(p/1000.);
	  if (lHit.energy()>0) p_mainParentEfrac_n->Fill(lHit.mainParentEfrac());
	  p_protonFrac_n->Fill(lHit.protonFrac());
	  p_neutronFrac_n->Fill(lHit.neutronFrac());
	  p_gFrac_n->Fill(lHit.gFrac());
	  p_eFrac_n->Fill(lHit.eFrac());
	  p_muFrac_n->Fill(lHit.muFrac());
	  p_hadFrac_n->Fill(lHit.hadFrac());
	  if (debug){
	    std::cout << " -- Found neutron! ";
	    lHit.Print(std::cout);
	  }
	  nNeutrons++;
	  break;
	}
      }

    }//loop on hits
      
  }//loop on entries
  
  std::cout << " -- Summary:" << std::endl
	    << " -- Found " << nNeutrons << " neutrons and " 
	    << nProtons << " protons and " 
	    << nPipm << " charged pions." << std::endl;    
  
  outputFile->Write();
  //outputFile->Close();
  
  return 0;
  
}//main
