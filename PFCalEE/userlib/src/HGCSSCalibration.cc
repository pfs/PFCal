#include "HGCSSCalibration.hh"
#include <sstream>
#include <iostream>
#include <cmath>

HGCSSCalibration::HGCSSCalibration(std::string filePath){

  vtx_x_ = 0;
  vtx_y_ = 0;
  vtx_z_ = 0;
  
  bool isScintOnly =  filePath.find("version22")!=filePath.npos;
  isHCALonly_ = filePath.find("version21")!=filePath.npos || isScintOnly;
  isCaliceHcal_ = filePath.find("version23")!=filePath.npos || filePath.find("version_23")!=filePath.npos;
  bool isECALHCAL = filePath.find("version20")!=filePath.npos;
  
  //set vertex pos for square section only
  if (filePath.find("model2")==filePath.npos){
    if (isCaliceHcal_) vtx_z_ = -1091.75;//mm
    else if (isScintOnly) vtx_z_ = -430.65;//mm
    else if (isHCALonly_) vtx_z_ = -824.01;//mm
    else if (isECALHCAL) vtx_z_ = -985.875;//mm
    else vtx_z_ = -161.865;//mm
  }

}


HGCSSCalibration::~HGCSSCalibration(){
}

double HGCSSCalibration::correctTime(const double & aTime,
				     const double & posx,
				     const double & posy,
				     const double & posz){
  double distance = sqrt(pow(posx-vtx_x_,2)+
			 pow(posy-vtx_y_,2)+
			 pow(posz-vtx_z_,2)
			 );
  double c = 299.792458;//3.e8*1000./1.e9;//in mm / ns...
  double cor = distance/c;
  // if (aTime>0 && cor > aTime) std::cout << " -- Problem ! Time correction is too large ";
  // if (aTime>0) std::cout << " -- hit time,x,y,z,cor = " 
  // 			 << aTime << " " << posx << " " << posy << " " << posz << " " 
  // 			 << cor << std::endl;
  double result = aTime-cor;
  if (result<0) result = 0;
  return result;
}

double HGCSSCalibration::MeVToMip(const unsigned layer, const bool absWeight) const{
  if (layer < theDetector().nLayers())
    return theDetector().subDetectorByLayer(layer).mipWeight
      *(absWeight?theDetector().subDetectorByLayer(layer).absWeight : 1.0);
  return 1;
}
