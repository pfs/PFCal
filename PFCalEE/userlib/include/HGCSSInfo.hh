#ifndef _hgcssinfo_hh_
#define _hgcssinfo_hh_
#include <iomanip>
#include <vector>
#include "Rtypes.h"
#include <sstream>
#include <map>
//#include "TH1F.h"

class HGCSSInfo{


public:
  HGCSSInfo()
  {
    version_ = -1;
    model_ = -1;
    cellsize_ = 2.5;//mm
    shape_ = 1;
  };
  
  ~HGCSSInfo(){};

  inline void version(const int & aVal){
    version_ = aVal;
  };
  
  inline int version() const{
    return version_;
  };

  inline void model(const int & aVal){
    model_ = aVal;
  };
  
  inline int model() const{
    return model_;
  };

  inline void shape(const unsigned & aVal){
    shape_ = aVal;
  };
  
  inline unsigned shape() const{
    return shape_;
  };

  inline void cellSize(const double & aVal){
    cellsize_ = aVal;
  };
  
  inline double cellSize() const{
    return cellsize_;
  };

  inline void calorSizeXY(const double & aVal){
    calorSizeXY_ = aVal;
  };
  
  inline double calorSizeXY() const{
    return calorSizeXY_;
  };

  inline void sensitiveZ(std::vector<double> &sensZ) {
    sensitiveZ_.clear();
    for(auto z : sensZ) sensitiveZ_.push_back(z);
  }
  inline const std::vector<double> &sensitiveZ() {
    return sensitiveZ_;
  }

  inline void etaBoundaryMin(std::vector<double> &etaMin) {
    etaBoundaryMin_.clear();
    for(auto z : etaMin) etaBoundaryMin_.push_back(z);
  }
  inline const std::vector<double> &etaBoundaryMin() {
    return etaBoundaryMin_;
  }

  inline void etaBoundaryMax(std::vector<double> &etaMax) {
    etaBoundaryMax_.clear();
    for(auto z : etaMax) etaBoundaryMax_.push_back(z);
  }
  inline const std::vector<double> &etaBoundaryMax() {
    return etaBoundaryMax_;
  }

private:

  int version_;
  int model_;
  double cellsize_;
  double calorSizeXY_;
  unsigned shape_;
  std::vector<double> sensitiveZ_,etaBoundaryMin_,etaBoundaryMax_;

  ClassDef(HGCSSInfo,3);
};


#endif
