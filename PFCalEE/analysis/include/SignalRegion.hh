#ifndef SignalRegion_h
#define SignalRegion_h

#include<string>
#include<iostream>
#include<fstream>
#include<sstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

#include "HGCSSEvent.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSPUenergy.hh"
#include "HGCSSGeometryConversion.hh"
#include "HGCSSCalibration.hh"
#include "PositionFit.hh"

#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"
#include "Math/Point3D.h"
#include "Math/Point3Dfwd.h"

class SignalRegion{

public:
    SignalRegion(const std::string inputFolder, 
                 const unsigned nLayers, 
                 const unsigned nevt,
                 const HGCSSGeometryConversion & geomConv,
                 const HGCSSPUenergy & puDensity,
		 const bool applyPuMixFix,
		 const unsigned versionNumber=12);

    ~SignalRegion();

  bool initialiseFitPositions();
  void initialise(TFile *outputFile,
		  const std::string outputDir);

  const FitResult & getAccurateFit(const unsigned ievt) const;

  ROOT::Math::XYZPoint getAccuratePos(const unsigned ievt, const unsigned iL) const;

  ROOT::Math::XYZPoint getAccuratePos(const FitResult & fit, const unsigned iL) const;

  Direction getAccurateDirection(const unsigned ievt) const;

  bool fillEnergies(const unsigned ievt,
		    const std::vector<HGCSSSamplingSection> & ssvec,
		    const std::vector<HGCSSSimHit> & simhitvec,
		    const std::vector<HGCSSRecoHit> & rechitvec,
		    const unsigned nPuVtx);

  bool fillEnergies(const unsigned ievt,
		    const HGCSSEvent & event,
		    const std::vector<HGCSSGenParticle> & genvec,
		    const std::vector<HGCSSSamplingSection> & ssvec,
		    const std::vector<HGCSSSimHit> & simhitvec,
		    const std::vector<HGCSSRecoHit> & rechitvec,
		    const unsigned nPuVtx);

  bool fillEnergies(const unsigned ievt,
		    const std::vector<HGCSSSamplingSection> & ssvec,
		    const std::vector<HGCSSSimHit> & simhitvec,
		    const std::vector<HGCSSRecoHit> & rechitvec,
		    const unsigned nPuVtx,
		    const FitResult & fit);

  bool fillEnergies(const unsigned ievt,
		    const std::vector<HGCSSSamplingSection> & ssvec,
		    const std::vector<HGCSSSimHit> & simhitvec,
		    const std::vector<HGCSSRecoHit> & rechitvec,
		    const unsigned nPuVtx,
		    const std::vector<ROOT::Math::XYZPoint> & eventPos);

  void finalise();
   
  void initialiseHistograms();
  void fillHistograms();
  
  inline void setOutputFile(TFile* outputFile){
    outputFile_ = outputFile;
    outputFile_->mkdir(outputDir_.c_str());
  };
  
  inline double getEtotalSR(const unsigned iSR, const bool subtractPU) const{
    double Etotal(0);
    if (iSR>=nSR_) return 0;
    for(unsigned iL(0);iL<nLayers_;iL++){
      if(subtractPU) Etotal += subtractedenergySR_[iL][iSR];
      else Etotal += energySR_[iL][iSR];
    }
    return Etotal;
  };
  
  inline double getSR(const unsigned iSR, const unsigned layer, const bool subtractPU) const{
    if(layer >= nLayers_) return 0;
    if (iSR>=nSR_) return 0;
    if(subtractPU) {return subtractedenergySR_[layer][iSR];}
    else {return energySR_[layer][iSR];}
  };

  inline double absweight(const unsigned layer) const{
    if(layer >= nLayers_) return 0;
    return absweight_[layer];
  };

  inline std::vector<std::vector<double> > getEnergyArray() const{
    return Exy_;
  };

private:
  
  unsigned nSR_;
  unsigned nevt_;
  std::string inputFolder_;
  unsigned nLayers_;    

  std::string outputDir_;
  TFile *outputFile_;
  TTree *outtree_;
  
  HGCSSGeometryConversion geomConv_;
  HGCSSPUenergy puDensity_;
  //HGCSSCalibration *mycalib_;
  
  bool fixForPuMixBug_;
  
  std::vector<FitResult> accurateFit_;

  unsigned nSkipped_;
  bool firstEvent_;

  //for tree
  std::vector<double> zPos_;
  std::vector<double> absweight_;
  unsigned evtIdx_;
  double totalE_;
  double wgttotalE_;
  double trueE_;
  double vtxX_;
  double vtxY_;
  double vtxZ_;
  double trueEta_;
  double truePhi_;
  std::vector<double> energyTot_;
  std::vector<std::vector<double> > energySR_;
  std::vector<std::vector<double> > subtractedenergySR_;
  std::vector<std::vector<double> > Exy_;

  TH1F *p_rawEtotal;
  TH1F *p_wgtEtotal;
  //std::vector<TH1F*> p_rawESR;
  std::vector<TH1F*> p_wgtESR;
  //std::vector<TH1F*> p_rawSubtractESR;
  std::vector<TH1F*> p_wgtSubtractESR;


};

#endif
