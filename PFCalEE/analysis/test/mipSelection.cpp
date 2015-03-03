#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#include <boost/algorithm/string.hpp>
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "boost/function.hpp"

// helpful tools
#include "KDTreeLinkerAlgoT.h"
#include <unordered_map>
#include <unordered_set>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
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

#include "TMath.h"
#include "TVectorD.h"
#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"
#include "Math/Point3D.h"
#include "Math/Point3Dfwd.h"
#include "TSystem.h"

#include "TVector3.h"

#include "utilities.h"

using boost::lexical_cast;
namespace po=boost::program_options;


struct MyHit{
  double x;
  double y;
  double z;
  double e;
  unsigned layer;
  bool hasSignal;
  ROOT::Math::XYZPoint position() const{
    return ROOT::Math::XYZPoint(x/10.,y/10.,z/10.);//in cm
  };
  double E() const {
    return e;
  };
};

void fillKDTree(KDTree & _hit_kdtree,
		std::vector<KDNode> & _hit_nodes,
		const std::vector<MyHit> & rechits){


  _hit_kdtree.clear();
  std::vector<bool> usable_rechits(rechits.size(),true);

  _hit_nodes.reserve(rechits.size());

  KDTreeCube kd_boundingregion = 
    fill_and_bound_kd_tree(rechits,usable_rechits,_hit_nodes);
  _hit_kdtree.build(_hit_nodes,kd_boundingregion);
  _hit_nodes.clear();
};

void findNeighbours(const double & cell_size_mm,
		    KDTree & _hit_kdtree,
		    const MyHit & current_cell, 
		    std::vector<KDNode> & found, 
		    const unsigned layerRange=1){
  const ROOT::Math::XYZPoint pos = current_cell.position();


  double cell_size = cell_size_mm/10.+0.01;//cm +0.01safety margin for double precision...
  unsigned layer = current_cell.layer;

  //std::cout << " -- x,y,z,l = " << pos.X() << " " << pos.Y() << " " << pos.Z() << " " << layer << std::endl;
  //std::cout << " -- kdtree size: " << _hit_kdtree.size() << std::endl;

  auto x_rh = minmax(pos.x()+cell_size,pos.x()-cell_size);
  auto y_rh = minmax(pos.y()+cell_size,pos.y()-cell_size);
  auto z_rh = minmax(zPos(layer-layerRange)/10.,zPos(layer+layerRange)/10.);//in cm!!!

  KDTreeCube hit_searchcube((float)x_rh.first,(float)x_rh.second,
			    (float)y_rh.first,(float)y_rh.second,
			    (float)z_rh.first,(float)z_rh.second);

  //std::cout << " -- min and max X = " << (float)x_rh.first << " " << (float)x_rh.second << std::endl
  //<< " -- min and max Y = " << (float)y_rh.first << " " <<(float)y_rh.second  << std::endl
  //<< " -- min and max Z = " << (float)z_rh.first << " " << (float)z_rh.second << std::endl;


  _hit_kdtree.search(hit_searchcube,found);

  //std::cout << " -- Found " << found.size() << " nodes." << std::endl;

};

bool matchingLayers(const MyHit & before2, 
		    const MyHit & before, 
		    const MyHit & central,
		    const MyHit & after,
		    const MyHit & after2){

  double precision = 0.1;

  if (fabs(before.x-central.x)<precision && 
      fabs(before.y-central.y)<precision &&
      fabs(after.x-central.x)<precision && 
      fabs(after.y-central.y)<precision &&
      fabs(before2.x-central.x)<precision && 
      fabs(before2.y-central.y)<precision &&
      fabs(after2.x-central.x)<precision && 
      fabs(after2.y-central.y)<precision) return true;

  /*if (fabs(before.x-central.x)<precision && before.y<central.y &&
      fabs(after.x-central.x)<precision && after.y>central.y) return true;
  if (fabs(before.x-central.x)<precision && before.y>central.y &&
      fabs(after.x-central.x)<precision && after.y<central.y) return true;
  if (before.x<central.x && fabs(before.y-central.y)<precision &&
      after.x>central.x && fabs(after.y-central.y)<precision) return true;
  if (before.x>central.x && fabs(before.y-central.y)<precision &&
      after.x<central.x && fabs(after.y-central.y)<precision) return true;

  if (before.x<central.x && before.y<central.y &&
      after.x>central.x && after.y>central.y) return true;
  if (before.x<central.x && before.y>central.y &&
      after.x>central.x && after.y<central.y) return true;
  if (before.x>central.x && before.y<central.y &&
      after.x<central.x && after.y>central.y) return true;
  if (before.x>central.x && before.y>central.y &&
  after.x<central.x && after.y<central.y) return true;*/
  
  return false;
};

int main(int argc, char** argv){//main  

  //std::string suffix = "_thresh0d8_notree";
  //std::string suffix = "_single";

  std::string cfg;
  double Ethresh;
  double EthreshMax;
  double EmaxCut;
  bool oneOnly;
  unsigned pNevts;
  std::string inFilePath;
  std::string outFilePath;
  unsigned debug;
  unsigned layerRange;
  unsigned start;

  po::options_description preconfig("Configuration"); 
  preconfig.add_options()("cfg,c",po::value<std::string>(&cfg)->required());
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(preconfig).allow_unregistered().run(), vm);
  po::notify(vm);
  po::options_description config("Configuration");
  config.add_options()
    ("pNevts,n", po::value<unsigned>(&pNevts)->default_value(0))
    ("inFilePath,i",   po::value<std::string>(&inFilePath)->required())
    ("outFilePath,o",  po::value<std::string>(&outFilePath)->required())
    ("Ethresh,e", po::value<double>(&Ethresh)->default_value(0.9))
    ("EthreshMax,m", po::value<double>(&EthreshMax)->default_value(5.0))
    ("EmaxCut", po::value<double>(&EmaxCut)->default_value(10000.))
    ("oneOnly",   po::value<bool>(&oneOnly)->default_value(true))
    ("debug,d",        po::value<unsigned>(&debug)->default_value(0))
    ("layerRange,l",   po::value<unsigned>(&layerRange)->default_value(1))
    ("start",        po::value<unsigned>(&start)->default_value(0))

    ;
  po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
  po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
  po::notify(vm);

 std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << inFilePath << std::endl
	    << " -- Output file path: " << outFilePath << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events." << std::endl;


  /////////////////////////////////////////////////////////////
  //hardcoded
  /////////////////////////////////////////////////////////////

  const double cell_size = 5;//mm

  const unsigned nEta = 1;
  //const unsigned nNoise = 10;
  const unsigned nNoise = 5;//10;

  //const double noise[nNoise] = {0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5};
  const double noise[nNoise] = {0.6,0.7,0.8,0.9,1.0};
  const double eta[nEta] = {2.85};//1.7,2.0,2.5};

  const double deta = 0.05;

  const unsigned nLayers = 30;

  /////////////////////////////////////////////////////////////
  //input
  /////////////////////////////////////////////////////////////

  //TFile *file = 0;
  //if (!testInputFile(inFilePath,file)) return 1;

  std::ostringstream inputbase;
  std::ostringstream input;
  TChain *tree = new TChain("MipStudy");
  unsigned nTrees = 1;  
  if (inFilePath.find(".root")!=inFilePath.npos) input << inFilePath;
  else {
    inputbase << inFilePath //<< "/MipStudy.root";
	      << "_run";
    nTrees = 2;
  }
  for (unsigned i(0); i<nTrees;++i){
    if (nTrees>1){
      input.str("");
      input << inputbase.str() << start+i
	    << ".root";
    }
    tree->AddFile(input.str().c_str());
    std::cout << "Adding MinBias file:" << input.str().c_str() << std::endl;
  }
  

  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  ///////// Output File // /////////////////////////
  //////////////////////////////////////////////////
  //////////////////////////////////////////////////


  std::ostringstream suffix;
  suffix << "Eta" << eta[0];
  suffix << "_thresh" << Ethresh ;
  suffix << "_" << EthreshMax;
  suffix << "_EmaxNeighbour" << EmaxCut;
  suffix << "_trk" << 1+2*layerRange << "layers";
  suffix << "_highnoise";
  if (cell_size > 6) suffix << "_1x1";
  if (oneOnly) suffix << "_onlyOne";
  if (nTrees>1) suffix << "_" << start << "_" << start+nTrees;

  std::ostringstream outPath;
  outPath << outFilePath << "/Histos" << suffix.str() << ".root";
  TFile *outputFile = TFile::Open(outPath.str().c_str(),"RECREATE");
  
  if (!outputFile) {
    std::cout << " -- Error, output file " << outPath.str() << " cannot be opened. Please create output directory. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
  }
  outputFile->cd();


  TH1F *hitSpectrum[nEta][nNoise][nLayers];
  TH1F *hitSpectrumSig[nEta][nNoise][nLayers];
  TH1F *p_nHits[nEta][nNoise][nLayers];

  TH1F *p_EmaxNeighbour[nEta][nNoise];

  for (unsigned ie(0);ie<nEta;++ie){
    for (unsigned in(0); in<nNoise;++in){
      std::ostringstream label;
      label.str("");
      label << "EmaxNeighbour_eta" << eta[ie]*10 << "_noise" << noise[in]*100;
      p_EmaxNeighbour[ie][in] = new TH1F(label.str().c_str(),";E_{max}^{neighbour} (Mips);hits",100,0,2);
      for (unsigned il(0); il<nLayers;++il){
	label.str("");
	label << "hitSpectrum_eta" << eta[ie]*10 << "_noise" << noise[in]*100 << "_layer" << il;
	hitSpectrum[ie][in][il] = new TH1F(label.str().c_str(),";E_{0.5#times0.5 cm^{2}} (mips);N_{cells}",100,0,5);
	label.str("");
	label << "hitSpectrumSig_eta" << eta[ie]*10 << "_noise" << noise[in]*100 << "_layer" << il;
	hitSpectrumSig[ie][in][il] = new TH1F(label.str().c_str(),";E_{0.5#times0.5 cm^{2}} (mips);N_{cells}",100,0,5);
	label.str("");
	label << "Nhits_eta" << eta[ie]*10 << "_noise" << noise[in]*100 << "_layer" << il;
	p_nHits[ie][in][il] = new TH1F(label.str().c_str(),";N_{0.5#times0.5 cm^{2}};N_{events}",500,0,5000);
      }
    }
  }

  //std::cout << " -- Getting tree" << std::endl;
  //file->cd();
  //TTree *tree = (TTree*)gDirectory->Get("MipStudy");

  std::vector<unsigned> * cellIdsX[nEta][nNoise];
  std::vector<unsigned> * cellIdsY[nEta][nNoise];
  std::vector<unsigned> * cellIdsZ[nEta][nNoise];
  std::vector<double> * energies[nEta][nNoise];
  std::vector<bool> * signal[nEta][nNoise];
  for (unsigned ie(0);ie<nEta;++ie){
    for (unsigned in(0); in<nNoise;++in){
      cellIdsX[ie][in] = 0;
      cellIdsY[ie][in] = 0;
      cellIdsZ[ie][in] = 0;
      energies[ie][in] = 0;
      signal[ie][in] = false;
      std::ostringstream label;
      label << "X_" << ie << "_" << in;
      tree->SetBranchAddress(label.str().c_str(),&cellIdsX[ie][in]);
      label.str("");
      label << "Y_" << ie << "_" << in;
      tree->SetBranchAddress(label.str().c_str(),&cellIdsY[ie][in]);
      label.str("");
      label << "Z_" << ie << "_" << in;
      tree->SetBranchAddress(label.str().c_str(),&cellIdsZ[ie][in]);
      label.str("");
      label << "E_" << ie << "_" << in;
      tree->SetBranchAddress(label.str().c_str(),&energies[ie][in]);
      label.str("");
      label << "signal_" << ie << "_" << in;
      tree->SetBranchAddress(label.str().c_str(),&signal[ie][in]);
    }
  }

  if (pNevts == 0) pNevts = tree->GetEntries();
  if (pNevts >tree->GetEntries()) pNevts = tree->GetEntries();

  std::cout << " -- Tree contains: " << tree->GetEntries() << " events." << std::endl;

  std::vector<KDNode> _hit_nodes;
  KDTree _hit_kdtree;

  for (unsigned ievt(0); ievt<pNevts; ++ievt){//loop on entries
    if (ievt%10 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
    tree->GetEntry(ievt);

    std::vector<MyHit> lHitVec[nEta][nNoise];

    for (unsigned ie(0);ie<nEta;++ie){
      for (unsigned in(0); in<nNoise;++in){
	for (unsigned ih(0);ih<(*cellIdsZ[ie][in]).size();++ih){
	  MyHit lHit;
	  lHit.e = (*energies[ie][in])[ih];
	  if (lHit.e==0) continue;
	  unsigned iL = (*cellIdsZ[ie][in])[ih];  
	  lHit.layer = iL;
	  int iX = (*cellIdsX[ie][in])[ih];
	  int iY = (*cellIdsY[ie][in])[ih];
	  lHit.x = -1695+ (iX-1)*cell_size + cell_size/2.;
	  lHit.y = -1695+ (iY-1)*cell_size + cell_size/2.;
	  lHit.z = zPos(iL);
	  lHit.hasSignal = (*signal[ie][in])[ih];
	  lHitVec[ie][in].push_back(lHit);
	}
      }
    }

    unsigned nHits[nEta][nNoise][nLayers];

    for (unsigned ie(0);ie<nEta;++ie){
      if (debug) std::cout << " - Eta = " << eta[ie] << std::endl;
      for (unsigned in(0); in<nNoise;++in){
	if (debug) std::cout << " --- Noise = " << noise[in] << std::endl;
	for (unsigned iL(0); iL<nLayers;++iL){
	  nHits[ie][in][iL] = 0;
	}

	fillKDTree(_hit_kdtree,
		   _hit_nodes,
		   lHitVec[ie][in]);
	
	if (debug) std::cout << " ----- Number of hits in vector : " << lHitVec[ie][in].size() << std::endl;
	//unsigned nMaxNeigh = 0;
	for (unsigned iH(0); iH<lHitVec[ie][in].size();++iH){
	  const MyHit& current = lHitVec[ie][in][iH];
	  unsigned iL = current.layer;

	  std::vector<KDNode> neighbours;
	  findNeighbours(cell_size,_hit_kdtree,current,neighbours,layerRange+1);
	  
	  if (debug>1) std::cout << " -- Number of closest neighbours found: " << neighbours.size() << std::endl; 
	  //if (neighbours.size()>nMaxNeigh) nMaxNeigh=neighbours.size();

	  bool passTrack = false;
	  double EmaxSC = 0;
	  std::vector<unsigned> idxBefore;
	  std::vector<unsigned> idxAfter;
	  std::vector<unsigned> idx2Before;
	  std::vector<unsigned> idx2After;
	  //clean neighbours
	  for( const KDNode& nbourpoint :neighbours ) {
	    unsigned index = nbourpoint.data;
	    if (index == iH) {
	      //neighbours.erase(neighbours.begin()+iH);
	      continue;
	    }
	    const MyHit& nbour = lHitVec[ie][in][index];
	    if ( nbour.layer != iL && (nbour.e <= Ethresh || nbour.e >= EthreshMax) ) {
	      //neighbours.erase(neighbours.begin()+iH);
	      continue;
	    }
	    if (nbour.layer == iL && (nbour.x==current.x || nbour.y==current.y )){
	      if (nbour.e>EmaxSC) EmaxSC = nbour.e;
	      //if (nbour.e>=EmaxCut) break;
	    }
	    else if ((iL>0 && nbour.layer == iL-1) || (iL==0 && nbour.layer == iL+2)) idxBefore.push_back(index);
	    else if ((iL<nLayers-1 && nbour.layer == iL+1) || (iL==nLayers-1 && nbour.layer == iL-2)) idxAfter.push_back(index);
	    else if (layerRange==2 && 
		     ((iL>1 && iL<nLayers-1 && nbour.layer == iL-2) || 
		      (iL==nLayers-1 && nbour.layer == iL-3) || 
		      (iL==1 && nbour.layer == iL+2) || 
		      (iL==0 && nbour.layer == iL+3))) idx2Before.push_back(index);
	    else if (layerRange==2 && 
		     ((iL>1 && iL<nLayers-2 && nbour.layer == iL+2) || 
		      (iL==0 && nbour.layer == iL+4) || 
		      (iL==1 && nbour.layer == iL+3) || 
		      (iL==nLayers-2 && nbour.layer == iL-3) || 
		      (iL==nLayers-1 && nbour.layer == iL-4))) idx2After.push_back(index);
	  }

	  const MyHit & ibef = idxBefore.size()>0?lHitVec[ie][in][idxBefore[0]]:current;
	  const MyHit & iaft = idxAfter.size()>0?lHitVec[ie][in][idxAfter[0]]:current;
	  const MyHit & i2bef = idx2Before.size()>0?lHitVec[ie][in][idx2Before[0]]:current;
	  const MyHit & i2aft = idx2After.size()>0?lHitVec[ie][in][idx2After[0]]:current;
	  if (!oneOnly) passTrack = idxBefore.size()>0 && idxAfter.size()>0  && ((layerRange==2 && idx2Before.size()>0 && idx2After.size()>0) || layerRange==1);
	  else passTrack = idxBefore.size()==1 && idxAfter.size()==1 && ((layerRange==2 && idx2Before.size()==1 && idx2After.size()==1) || layerRange==1) && matchingLayers(i2bef,ibef,current,iaft,i2aft);
	  if (!passTrack) continue;

	  p_EmaxNeighbour[ie][in]->Fill(EmaxSC);
	  if (EmaxSC>=EmaxCut) continue;

	  hitSpectrum[ie][in][iL]->Fill(current.e);
	  if (current.hasSignal) hitSpectrumSig[ie][in][iL]->Fill(current.e);
	  nHits[ie][in][iL]++;
	  
	}//loop on hits
	//std::cout << " -- Max N neighbours = " << nMaxNeigh << std::endl; 
	for (unsigned iL(0); iL<nLayers;++iL){
	  p_nHits[ie][in][iL]->Fill(nHits[ie][in][iL]);
	}//loop on layers
	
      }//loop on noise
    }//loop on eta

  }//loop on tree

  outputFile->cd();
  outputFile->Write();
  outputFile->Close();

  return 0;

}//main
