#include "HGCSSCalibration.hh"
#include <sstream>
#include <iostream>
#include <cmath>

HGCSSCalibration::HGCSSCalibration(std::string filePath,
				   const bool bypassR,
				   const unsigned nSi){

  vtx_x_ = 0;
  vtx_y_ = 0;
  vtx_z_ = 0;
  bypassRadius_ = bypassR;
  nSiLayers_ = nSi;
}


HGCSSCalibration::~HGCSSCalibration(){
}

double HGCSSCalibration::addTimeOfFlight(const double & aTime,
					 const double & posx,
					 const double & posy,
					 const double & posz){
  return addTimeOfFlight(aTime,posx,posy,posz,vtx_x_,vtx_y_,vtx_z_);
}

double HGCSSCalibration::addTimeOfFlight(const double & aTime,
					 const double & posx,
					 const double & posy,
					 const double & posz,
					 const double & vtxx,
					 const double & vtxy,
					 const double & vtxz){
  double distance = sqrt(pow(posx-vtxx-vtx_x_,2)+
			 pow(posy-vtxy-vtx_y_,2)+
			 pow(posz-vtxz-vtx_z_,2)
			 );
  double c = 299.792458;//3.e8*1000./1.e9;//in mm / ns...
  double cor = distance/c;
  // if (aTime>0 && cor > aTime) std::cout << " -- Problem ! Time correction is too large ";
  // if (aTime>0) std::cout << " -- hit time,x,y,z,cor = " 
  // 			 << aTime << " " << posx << " " << posy << " " << posz << " " 
  // 			 << cor << std::endl;
  double result = aTime+cor;
  return result;
}

double HGCSSCalibration::correctTime(const double & aTime,
				     const double & posx,
				     const double & posy,
				     const double & posz){
  return correctTime(aTime,posx,posy,posz,vtx_x_,vtx_y_,vtx_z_);
}

double HGCSSCalibration::correctTime(const double & aTime,
				     const double & posx,
				     const double & posy,
				     const double & posz,
				     const double & vtxx,
				     const double & vtxy,
				     const double & vtxz){
  double distance = sqrt(pow(posx-vtxx,2)+
			 pow(posy-vtxy,2)+
			 pow(posz-vtxz,2)
			 );
  double c = 299.792458;//3.e8*1000./1.e9;//in mm / ns...
  double cor = distance/c;
  // if (aTime>0 && cor > aTime) std::cout << " -- Problem ! Time correction is too large ";
  // if (aTime>0) std::cout << " -- hit time,x,y,z,cor = " 
  // 			 << aTime << " " << posx << " " << posy << " " << posz << " " 
  // 			 << cor << std::endl;
  double result = aTime-cor;
  //if (result<0) result = 0;
  return result;
}

double HGCSSCalibration::MeVToMip(const unsigned layer, const bool absWeight) const{
  if (layer < theDetector().nLayers())
    return theDetector().subDetectorByLayer(layer).mipWeight
      *(absWeight?theDetector().subDetectorByLayer(layer).absWeight : 1.0);
  return 1;
}

/*double HGCSSCalibration::MeVToMip(const unsigned layer, const double aEta, const bool absWeight) const{
  double res = 1;
  if (layer < theDetector().nLayers())
    res = theDetector().subDetectorByLayer(layer).mipWeight
      *(absWeight?theDetector().subDetectorByLayer(layer).absWeight : 1.0);

  if (aEta<=1.75) return res*2./3.;//300um
  else if (aEta > 2.15) return res*2.;//100um
  return res;

  }*/

double HGCSSCalibration::MeVToMip(const unsigned layer, const double aRadius, const bool absWeight) const{
  double res = 1;
  if (layer < theDetector().nLayers())
    res = theDetector().subDetectorByLayer(layer).mipWeight
      *(absWeight?theDetector().subDetectorByLayer(layer).absWeight : 1.0);

  if (theDetector().subDetectorByLayer(layer).isSi == false) return res;

  if (bypassRadius_) return res*2./nSiLayers_;

  double r1 = 1200;
  double r2 = theDetector().subDetectorByLayer(layer).radiusLim;
  if (theDetector().subDetectorByLayer(layer).type == DetectorEnum::FHCAL) {
    r1 = 1000;
  }
  if (aRadius>r1) return 1./0.0813;//res*2./3.;//300um
  else if (aRadius < r2) return 1./0.0246;//res*2.;//100um
  return 1./0.0524;//res;
}
