#ifndef _hgcsscluster_hh_
#define _hgcsscluster_hh_

#include <iomanip>
#include <vector>
#include <cmath>
#include "Rtypes.h"
#include <sstream>
#include "TMath.h"

#include "HGCSSRecoHit.hh"
#include "Math/VectorUtil.h"
#include "Math/Point3D.h"
#include "Math/Point3Dfwd.h"
#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "TVector3.h"
#include "Math/GenVector/VectorUtil.h"
//#include "Math/GenVector/BitReproducible.h"

class HGCSSCluster{

public:
  HGCSSCluster():
    energy_(0),
    pos_(0,0,0),
    dir_(0,0,0),
    layer_(0),
    seedPos_(0,0,0),
    seedE_(0),
    width_(0)
  {
  };

  HGCSSCluster(const HGCSSRecoHit & aRecHit);

  ~HGCSSCluster(){};

  //double eta() const;
  //double theta() const;
  //double phi() const;

  inline double getSeedE() const{
    return seedE_;
  };

  double getSeedEta() const;
  double getSeedPhi() const;

  inline ROOT::Math::XYZPoint position() const{
    return pos_;
  };

  inline ROOT::Math::XYZVector direction() const{
    return dir_;
  };

   inline ROOT::Math::XYZPoint seedPosition() const{
    return seedPos_;
  };
  
  inline double energy() const {
    return energy_;
  };

  inline double pt() const {
    return energy_/cosh(dir_.eta());
  };

  inline double px() const {
    return pt()*cos(dir_.phi());
  };

  inline double py() const {
    return pt()*sin(dir_.phi());
  };

  inline double pz() const {
    return pt()*sinh(dir_.eta());
  };

  inline void setEnergy(const double & energy) {
    energy_ = energy;
  };

  inline void setSeedEnergy(const double & energy) {
    seedE_ = energy;
  };

  inline void setPosition(const ROOT::Math::XYZPoint & pos){
    pos_ = pos;
  };

  inline void setDirection(const ROOT::Math::XYZVector & dir){
    dir_ = dir;
  };

  inline unsigned layer() const {
    return layer_;
  };

  inline unsigned width() const{
    return width_;
  };

  inline void setLayer(const unsigned & layer){
    layer_ = layer;
  };

  inline void setSeed(const ROOT::Math::XYZPoint & pos){
    seedPos_ = pos;
  };

  void addRecHitFraction(std::pair<HGCSSRecoHit*,double> aHit);
  const std::map<HGCSSRecoHit*,double> & recHitFractions() const;

  inline unsigned nRecHits() const{
    return recHitMap_.size();
  }

  void calculatePosition();

  void calculateDirection();
  void calculateDirection(const std::vector<unsigned> & lToRemove);


  void Print(std::ostream & aOs) const;

private:

  double energy_;
  ROOT::Math::XYZPoint pos_;
  ROOT::Math::XYZVector dir_;
  unsigned layer_;
  ROOT::Math::XYZPoint seedPos_;
  double seedE_;
  unsigned width_;
  std::map<HGCSSRecoHit*,double> recHitMap_;

  ClassDef(HGCSSCluster,1);

};


typedef std::vector<HGCSSCluster> HGCSSClusterVec;



#endif
