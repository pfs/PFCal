#ifndef PositionFit_hh
#define PositionFit_hh

#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"

#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSCluster.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSPUenergy.hh"
#include "HGCSSGeometryConversion.hh"

#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"

#include "TVector3.h"

struct FitResult{
  double pos_x;
  double pos_y;
  double tanangle_x;
  double tanangle_y;
  bool found;
  FitResult():pos_x(0),pos_y(0),tanangle_x(0),tanangle_y(0),found(false)
  {};
};

struct Direction{
  double tanangle_x;
  double tanangle_y;
  double eta() const{
    return asinh(1.0/sqrt(tanangle_x*tanangle_x+tanangle_y*tanangle_y));
  };
  double phi() const{
    return atan2(tanangle_y,tanangle_x);
  };
  TVector3 dir() const{
    TVector3 v;
    v.SetPtEtaPhi(1./cosh(eta()),eta(),phi());
    return v;
  };
  
  Direction(){};
  Direction(double tanx,double tany){
    tanangle_x = tanx;
    tanangle_y = tany;
  };
  double GetX(double posz,double vtx_x, double vtxz)const{
    return vtx_x+tanangle_x*(posz-vtxz);
  };
  double GetY(double posz,double vtx_y, double vtxz)const{
    return vtx_y+tanangle_y*(posz-vtxz);
  };
  ROOT::Math::XYZPoint GetPosFF(const ROOT::Math::XYZPoint & vtx)const{
    double zposFF = 3173.9;
    return ROOT::Math::XYZPoint(GetX(zposFF,vtx.x(),vtx.z()),GetY(zposFF,vtx.y(),vtx.z()),zposFF);
  };
  void Print() const{
    std::cout << "tanangles = " << tanangle_x << " " << tanangle_y << " pt=" << dir().Pt() << " eta=" << eta() << " " << dir().Eta() << " phi=" << phi() << " " << dir().Phi() << " p=" << dir().Mag() << std::endl;
  };
};


class PositionFit{

public:

  PositionFit(const unsigned nSR,
	      const double & residualMax, 
	      const unsigned nLayers, 
	      const unsigned nSiLayers,
	      const bool applyPuMixFix,
	      const unsigned debug=0,
	      const bool doMatrix=true);

  ~PositionFit(){};

  double getW0(const unsigned layer);

  double DeltaPhi(const double & phi1, const double & phi2);

  std::pair<unsigned, std::pair<double,double> > findMajorityValue(std::vector<std::pair<double,double> > & values) const;

  void initialise(TFile* outputFile,
		  const std::string outputDir, 
		  const std::string outFolder, 
		  const HGCSSGeometryConversion & geomConv, 
		  const HGCSSPUenergy & puDensity);

  inline void setMatrixFolder(const std::string outFolder){
    matrixFolder_ = outFolder;
  };

  void initialiseClusterHistograms();
  void initialisePositionHistograms();
  void initialiseFitHistograms();

  bool getZpositions(const unsigned versionNumber);
  void getZpositions(const unsigned versionNumber,
		     TTree *aSimTree,
		     const unsigned nEvts);

  void getInitialPositions(TTree *simTree, 
			   TTree *recoTree,
			   const unsigned nEvts,
			   const unsigned G4TrackID=1);

  bool getInitialPosition(const unsigned ievt,
			  const unsigned nVtx, 
			  std::vector<HGCSSRecoHit> *rechitvec,
			  unsigned & nTooFar,
			  std::vector<std::vector<double> > & Exy);

  bool getGlobalMaximum(const unsigned ievt, 
			const unsigned nVtx, 
			std::vector<HGCSSRecoHit> *rechitvec, 
			const ROOT::Math::XYZVector & truthPos0, 
			double & phimax,double & etamax);

  void findSeeds(std::vector<HGCSSRecoHit> *rechitvec,
		 std::vector<bool> & seedable);

  unsigned getClusters(std::vector<HGCSSRecoHit> *rechitvec,
		       HGCSSClusterVec & output);

  bool setTruthInfo(std::vector<HGCSSGenParticle> *genvec, const int G4TrackID);

  //bool getTruthPosition(std::vector<HGCSSGenParticle> *genvec,std::vector<ROOT::Math::XYPoint> & truthPos, const int trackID=1);

  void getMaximumCellFromGeom(const double & phimax,const double & etamax,std::vector<double> & xmax,std::vector<double> & ymax);

  void getMaximumCell(std::vector<HGCSSRecoHit> *rechitvec,const double & phimax,const double & etamax,std::vector<double> & xmax,std::vector<double> & ymax);

  void getEnergyWeightedPosition(std::vector<HGCSSRecoHit> *rechitvec,
				 const unsigned nPU, 
				 const std::vector<double> & xmax,
				 const std::vector<double> & ymax,
				 std::vector<ROOT::Math::XYPoint> & recoPos,
				 std::vector<double> & recoE,
				 std::vector<unsigned> & nHits,
				 std::vector<double> & puE,
				 std::vector<std::vector<double> > & Exy,
				 const bool puSubtracted=true);

  void getPuContribution(std::vector<HGCSSRecoHit> *rechitvec, const std::vector<double> & xmax,const std::vector<double> & ymax,std::vector<double> & puE);

  void fillErrorMatrix(const std::vector<ROOT::Math::XYPoint> & recoPos, const std::vector<unsigned> & nHits);

  void finaliseErrorMatrix(const bool doX);
  void finaliseErrorMatrix();
  bool fillMatrixFromFile(const bool doX, const bool old);
  bool fillMatrixFromFile(const bool old=false);
  void fillCorrelationMatrix();

  bool getPositionFromFile(const unsigned ievt,
			   std::vector<unsigned> & layerId,
			   std::vector<double> & posx,
			   std::vector<double> & posy,
			   std::vector<double> & posz,
			   std::vector<double> & posxtruth,
			   std::vector<double> & posytruth,
			   std::vector<double> & E,
			   bool cutOutliers=false,
			   bool print=true);

  bool initialiseLeastSquareFit();
  //return 1 if no input file or <3 layers
  //return 2 if chi2/ndf>chi2ndfmax_
  //return 0 if success
  unsigned performLeastSquareFit(const unsigned ievt,
				 FitResult & fit);
  void finaliseFit();

  //return 1 if no input file or <3 layers
  //return 2 if chi2/ndf>chi2ndfmax_
  //return 0 if success
  unsigned fitEvent(const unsigned ievt,
		    FitResult & fit,
		    const bool cutOutliers=false);

  unsigned GetIndependentFitResult(const TMatrixDSym & e,
				   const unsigned nL,
				   //const unsigned rt,
				   const unsigned ievt,
				   const TVectorD & z,
				   const TVectorD & x,
				   const TVectorD & y,
				   const TVectorD & xy,
				   TMatrixD & w,
				   TVectorD & p,
				   double & chiSq
				   );

  unsigned GetSimultaneousFitResult(const TMatrixDSym & e,
				const unsigned nL,
				//const unsigned rt,
				const unsigned ievt,
				const TVectorD & z,
				const TVectorD & xy,
				TVectorD & simPars,
				TVectorD & delta,
				double & chiSq
				);

  inline void setOutputFile(TFile * outputFile){
    outputFile_ = outputFile;
    outputFile_->mkdir(outputDir_.c_str());
  };

  inline unsigned nSR() const{
    return nSR_;
  };

  inline double residualMax() const{
    return residualMax_;
  };

  inline void setRecoDir(double tanx,double tany){
    recoDir_ = Direction(tanx,tany);
  };

  inline void setTruthDir(double tanx,double tany){
    truthDir_ = Direction(tanx,tany);
  };

  inline const Direction & recoDir() const{
    return recoDir_;
  };

  inline const Direction & truthDir() const{
    return truthDir_;
  };

  inline const ROOT::Math::XYZPoint & truthVtx() const{
    return truthVtx_;
  };

  inline double truthE() const{
    return truthE_;
  };

  inline ROOT::Math::XYPoint truthPos(const unsigned iL) const{
    return ROOT::Math::XYPoint(truthDir_.GetX(avgZ_[iL],truthVtx_.x(),truthVtx_.z()),
			       truthDir_.GetY(avgZ_[iL],truthVtx_.y(),truthVtx_.z()));

  };

private:
  PositionFit(){};

  unsigned nSR_;
  double residualMax_;
  double chi2ndfmax_;
  double precision_;
  double seedMipThreshold_;
  double maxdR_;
  unsigned nLayers_;
  unsigned nSiLayers_;
  double cellSize_;
  unsigned debug_;
  bool useMeanPU_;
  bool fixForPuMixBug_;
  bool doMatrix_;
  bool saveEtree_;
  bool doLogWeight_;

  HGCSSGeometryConversion geomConv_;
  HGCSSPUenergy puDensity_;

  std::vector<double> avgZ_;

  //initialisation for error matrix
  std::vector<double> mean_;//sum residuals for x and y
  std::vector<std::vector<double> > sigma_;//sum square
  std::vector<unsigned> nL_mean_;//number of valid layers
  std::vector<std::vector<unsigned> > nL_sigma_;
  TMatrixD matrix_;
  TMatrixD corrMatrix_;

  Direction recoDir_;
  Direction truthDir_;
  ROOT::Math::XYZPoint truthVtx_;
  double truthE_;

  unsigned nInvalidFits_;
  unsigned nFailedFitsAfterCut_;
  unsigned countFailedSimFits_;
  std::ofstream fout_;

  //path for saving data files
  std::string outFolder_;
  std::string matrixFolder_;
  std::string outputDir_;
  TFile *outputFile_;

  //cluster histos
  TH1F *p_nClusters;

  TH1F *p_clusnHits_all;
  TH1F *p_seedEoverE_all;
  TH1F *p_clusLayer_all;
  TH1F *p_clusWidth_all;
  TH1F *p_seeddeta_all;
  TH1F *p_seeddphi_all;

  TH1F *p_clusnHits_sel;
  TH1F *p_seedEoverE_sel;
  TH1F *p_clusLayer_sel;
  TH1F *p_clusWidth_sel;
  TH1F *p_seeddeta_sel;
  TH1F *p_seeddphi_sel;

  TH1F *p_mindRtruth;

  //position histos
  TH1F *p_nGenParticles;

  std::vector<TH2F *> p_genxy;
  std::vector<TH2F *> p_recoxy;
  std::vector<TH1F *> p_dRmin;

  //TH1F *p_numberOfMaxTried;
  //TH1F *p_dRMaxTruth;
  TH2F *p_hitEventPuContrib;
  TH2F *p_hitMeanPuContrib;
  TH1F *p_diffPuContrib;

  TH1F *p_residuals_x;
  TH1F *p_residuals_y;
  TH2F *p_etavsphi;
  TH2F *p_etavsphi_max;
  TH2F *p_etavsphi_truth;

  TH2F *p_yvsx_max;
  TH2F *p_yvsx_truth;

  //fit histos
  TH1F *p_iterations;
  TH1F *p_nLayersFit;
  TH2F *p_recoXvsLayer;
  TH2F *p_recoYvsLayer;
  TH2F *p_recoZvsLayer;
  TH2F *p_truthXvsLayer;
  TH2F *p_truthYvsLayer;
  TH2F *p_fitXvsLayer;
  TH2F *p_fitYvsLayer;
  TH2D *p_errorMatrix_xy;
  TH2D *p_corrMatrix_xy;
  //TH2D *p_errorMatrix_y;
  //TH2D *p_corrMatrix_y;
  TH1F *p_chi2[2];
  TH1F *p_chi2overNDF[2];
  TH1F *p_impactX0[2];
  TH1F *p_impactY0[2];
  TH1F *p_impactXFF[2];
  TH1F *p_impactYFF[2];
  TH1F *p_impactX14[2];
  TH1F *p_impactY14[2];
  TH1F *p_impactZx[2];
  TH1F *p_impactZy[2];
  TH1F *p_tanAngleX[2];
  TH1F *p_tanAngleY[2];
  TH1F *p_positionReso;
  TH1F *p_angularReso;
  TH1F *p_impactX0_residual;
  TH1F *p_impactY0_residual;
  TH1F *p_impactXFF_residual;
  TH1F *p_impactYFF_residual;
  TH1F *p_impactX14_residual;
  TH1F *p_impactY14_residual;
  TH1F *p_impactZx_residual;
  TH1F *p_impactZy_residual;
  TH1F *p_tanAngleX_residual;
  TH1F *p_tanAngleY_residual;
  TH1F *p_angleX_residual;
  TH1F *p_angleY_residual;
  TH1F *p_eta_reco;
  TH1F *p_phi_reco;
  TH1F *p_eta_truth;
  TH1F *p_phi_truth;
  TH1F *p_eta_residual;
  TH1F *p_phi_residual;

};//class

#endif
