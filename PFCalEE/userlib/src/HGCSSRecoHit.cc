#include "HGCSSRecoHit.hh"

#include <iomanip>
#include <cmath>
#include <stdlib.h>

HGCSSRecoHit::HGCSSRecoHit(const HGCSSSimHit & aSimHit, const unsigned granularity){
  energy_ = aSimHit.energy();
  adcCounts_ = 0;
  zpos_ = aSimHit.get_z();


  layer_ = aSimHit.layer();
  noiseFrac_ = 0;

  xpos_ = aSimHit.get_x();
  ypos_ = aSimHit.get_y();
  //cellid encoding:
  //bool x_side = x>0 ? true : false;
  //bool y_side = y>0 ? true : false;
  //unsigned x_cell = static_cast<unsigned>(fabs(x)/(cellSize*granularity));
  //unsigned y_cell = static_cast<unsigned>(fabs(y)/(cellSize*granularity));
  //encodeCellId(x_side,y_side,x_cell,y_cell,granularity);

}


double HGCSSRecoHit::theta() const {
  return 2*atan(exp(-eta()));
}

double HGCSSRecoHit::eta() const {
  double theta = acos(fabs(zpos_)/sqrt(zpos_*zpos_+xpos_*xpos_+ypos_*ypos_));
  double leta = -log(tan(theta/2.));
  if (zpos_>0) return leta;
  else return -leta;
}

double HGCSSRecoHit::phi() const {
  double x = get_x();
  double y = get_y();
  if (x==0) return 0;
  if (x>0) return atan(y/x);
  else if (y>0) return TMath::Pi()+atan(y/x);
  else return -TMath::Pi()+atan(y/x);
}

// void HGCSSRecoHit::encodeCellId(const bool x_side,const bool y_side,const unsigned x_cell,const unsigned y_cell, const unsigned granularity){
//   cellid_ = 
//     x_side | ((x_cell & 0xFFF)<<1) |
//     (y_side<<13) | ((y_cell & 0xFFF)<<14) |
//     ((granularity & 0x3F) <<26) ;

//   // std::cout << " Cross-check of encoding: cellid=" << cellid_ << std::endl
//   // 	    << " x_side " << x_side << " " << get_x_side() << std::endl
//   // 	    << " y_side " << y_side << " " << get_y_side() << std::endl
//   // 	    << " x_cell " << x_cell << " " << get_x_cell() << std::endl
//   // 	    << " y_cell " << y_cell << " " << get_y_cell() << std::endl
//   //        << " granularity " << granularity << " " << getGranularity() << std::endl
//   //   ;
// }

void HGCSSRecoHit::Add(const HGCSSSimHit & aSimHit){
  energy_ += aSimHit.energy();
}

void HGCSSRecoHit::Print(std::ostream & aOs) const{
  aOs << "====================================" << std::endl
      << " = Layer " << layer_ //<< " cellid " << cellid_ 
      << std::endl
      << " = Energy " << energy_ << " noiseFrac " << noiseFrac_ << std::endl
      << " = Digi E " << adcCounts_ << " adcCounts." << std::endl
      << "====================================" << std::endl;

}
