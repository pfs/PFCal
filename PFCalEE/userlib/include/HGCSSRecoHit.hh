#ifndef _hgcssrecohit_hh_
#define _hgcssrecohit_hh_

#include <iomanip>
#include <vector>
#include <cmath>
#include "Rtypes.h"
#include <sstream>
#include "TMath.h"

#include "HGCSSSimHit.hh"

class HGCSSRecoHit{

public:
  HGCSSRecoHit():
    energy_(0),
    adcCounts_(0),
    xpos_(0),
    ypos_(0),
    zpos_(0),
    layer_(0),
    //cellid_(0),
    noiseFrac_(0)
  {};

  HGCSSRecoHit(const HGCSSSimHit & aSimHit, const unsigned granularity=1);

  ~HGCSSRecoHit(){};

  inline double energy() const {
    return energy_;
  };

  double eta() const;
  double theta() const;
  double phi() const;

  inline double E() const {
    return energy_;
  };

  inline double pt() const {
    return energy_/cosh(eta());
  };

  inline double px() const {
    return pt()*cos(phi());
  };

  inline double py() const {
    return pt()*sin(phi());
  };

  inline double pz() const {
    return pt()*sinh(eta());
  };

  inline void energy(const double & energy) {
    energy_ = energy;
  };

  inline void x(const double & pos){
    xpos_ = pos;
  };

  inline void y(const double & pos){
    ypos_ = pos;
  };

  inline void z(const double & pos){
    zpos_ = pos;
  };

  inline unsigned adcCounts() const {
    return adcCounts_;
  };

  inline void adcCounts(const unsigned & adcCounts){
    adcCounts_ = adcCounts;
  };

  inline unsigned layer() const {
    return layer_;
  };

  inline void layer(const unsigned & layer){
    layer_ = layer;
  };

  // inline unsigned cellid() const {
  //   return cellid_;
  // };

  // inline void cellid(const unsigned & id){
  //   cellid_ = id;
  // };

  inline double noiseFraction() const {
    return noiseFrac_;
  };

  inline void noiseFraction(const double & aFrac){
    noiseFrac_ = aFrac;
  };

  void Add(const HGCSSSimHit & aSimHit);
  
  // void encodeCellId(const bool x_side,const bool y_side,const unsigned x_cell,const unsigned y_cell, const unsigned granularity);

  // inline bool get_x_side() const{
  //   return cellid_ & 0x0001;
  // };

  // inline bool get_y_side() const{
  //   return (cellid_ & 0x2000) >> 13;
  // };

  // inline unsigned get_x_cell() const{
  //   return (cellid_ & 0x1FFE) >> 1;
  // };

  // inline unsigned get_y_cell() const{
  //   return (cellid_ & 0x03FFC000) >> 14;
  // };

  inline double get_x() const{
    return xpos_;
    //float sign = get_x_side() ? 1. : -1. ;
    //if (sign > 0)
    //  return get_x_cell()*sign*cellSize*getGranularity()+cellSize*getGranularity()/2;
    //else return get_x_cell()*sign*cellSize*getGranularity()-cellSize*getGranularity()/2;
  };

  inline double get_y() const{
    return ypos_;
    //float sign = get_y_side() ? 1. : -1. ;
    //if (sign > 0)
    //  return get_y_cell()*sign*cellSize*getGranularity()+cellSize*getGranularity()/2;
    //else return get_y_cell()*sign*cellSize*getGranularity()-cellSize*getGranularity()/2;
  };

  inline double get_z() const {
    return zpos_;
  };

  // inline unsigned getGranularity() const{
  //   return (cellid_ & 0xFC000000) >> 26;
  // };

  void Print(std::ostream & aOs) const;

private:

  double energy_;
  unsigned adcCounts_;
  double xpos_;
  double ypos_;
  double zpos_;
  unsigned layer_;
  //unsigned cellid_;
  double noiseFrac_;

  ClassDef(HGCSSRecoHit,1);

};


typedef std::vector<HGCSSRecoHit> HGCSSRecoHitVec;



#endif
