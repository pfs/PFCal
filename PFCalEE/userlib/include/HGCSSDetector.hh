#ifndef HGCSSDetector_h
#define HGCSSDetector_h

#include <string>
#include <vector>
#include <map>
#include "TH2D.h"

enum DetectorEnum {
  FECAL,
  MECAL,
  BECAL,
  FHCAL,
  BHCAL1,
  BHCAL2
};

class HGCSSSubDetector {

public:
  HGCSSSubDetector():
    type(DetectorEnum::FECAL),
    name(""),
    layerIdMin(0),
    layerIdMax(0),
    mipWeight(1),
    absWeight(1),
    gevWeight(1),
    gevOffset(0),
    isSi(false),
    isScint(false),
    radiusLim(0)
  {};
  ~HGCSSSubDetector(){};

  DetectorEnum type;
  std::string name;
  unsigned layerIdMin;
  unsigned layerIdMax;
  double mipWeight;
  double absWeight;
  double gevWeight;
  double gevOffset;
  bool isSi;
  bool isScint;
  double radiusLim;

  inline unsigned nLayers() const{
    return (layerIdMax-layerIdMin);
  };

private:

};

class HGCSSDetector {

public:
  friend HGCSSDetector & theDetector();

  inline void initialiseIndices(const unsigned versionNumber){
    
    indices_.clear();
    indices_.resize(7,0);
    //fill layer indices
    if (versionNumber==22){
      indices_[4] = 0;
      indices_[5] = 10;
      indices_[6] = 10;
    }
    else if (versionNumber==32) {
      indices_[4] = 0;
      indices_[5] = 12;
      indices_[6] = 12;
    }
    else if (versionNumber==23) {
      indices_[3] = 0;
      indices_[4] = 38;
      indices_[5] = 47;
      indices_[6] = 54;
    }
    else if (versionNumber==21) {
      indices_[3] = 0;
      indices_[4] = 24;
      indices_[5] = 34;
      indices_[6] = 34;
    }
    else if (versionNumber==31) {
      indices_[3] = 0;
      indices_[4] = 12;
      indices_[5] = 24;
      indices_[6] = 24;
    }
    else if (versionNumber < 20){
      indices_[0] = 0;
      indices_[1] = versionNumber==8?11:10;
      indices_[2] = versionNumber==8?21:20;
      indices_[3] = versionNumber==8?31:30;
      indices_[4] = indices_[3];
      indices_[5] = indices_[3];
      indices_[6] = indices_[3];
    }
    else if (versionNumber == 30){
      indices_[0] = 0;
      indices_[1] = 10;
      indices_[2] = 20;
      indices_[3] = 29; // For Si pre-shower
      indices_[4] = indices_[3];
      indices_[5] = indices_[3];
      indices_[6] = indices_[3];
    }
    else if (versionNumber == 33){
      indices_[0] = 0;
      indices_[1] = 10;
      indices_[2] = 20;
      indices_[3] = 28;
      indices_[4] = 40;
      indices_[5] = 52;
      indices_[6] = 52;
    }
    else {
      indices_[0] = 0;
      indices_[1] = versionNumber==24?11:10;
      indices_[2] = versionNumber==24?21:20;
      indices_[3] = versionNumber==24?31:30;
      indices_[4] = versionNumber==24?55:54;
      indices_[5] = versionNumber==24?65:66;
      indices_[6] = versionNumber==24?65:66;
    }
    
  };

  void buildDetector(const unsigned versionNumber,
		     bool concept=true,
		     bool isCaliceHcal=false,
		     bool bypassR=false);

  const HGCSSSubDetector & subDetectorByLayer(const unsigned aLayer);

  unsigned getSection(const unsigned aLayer) const;
  inline unsigned section(const DetectorEnum adet){
    if (enumMap_.find(adet) != enumMap_.end())
      return enumMap_[adet];
    return nSections();
  };

  void addSubdetector(const HGCSSSubDetector & adet);
  
  void finishInitialisation();

  inline unsigned nLayers(const unsigned aSection) const{
    return subdets_[aSection].nLayers();
  };

  inline unsigned nLayers(DetectorEnum adet){
    return subdets_[enumMap_[adet]].nLayers();
  };

  const HGCSSSubDetector & subDetectorByEnum(DetectorEnum adet);
  inline const HGCSSSubDetector & subDetectorBySection(const unsigned aSection) const{
    return subdets_[aSection];
  };

  inline unsigned nLayers() const{
    return nLayers_;
  };

  inline unsigned nSections() const{
    return nSections_;
  };

  inline DetectorEnum detType(const unsigned aSection) const{
    return subdets_[aSection].type;
  };
  
  inline DetectorEnum detTypeLayer(const unsigned aLayer) const{
    return subdets_[getSection(aLayer)].type;
  };
  
  inline std::string detName(const unsigned aSection) const{
    return subdets_[aSection].name;
  };


  void reset();

  void printDetector(std::ostream & aOs) const ;

private:
  HGCSSDetector(){
    bypassRadius_ = false;
  };

  ~HGCSSDetector(){
    reset();
  };
  
  std::vector<HGCSSSubDetector> subdets_;
  std::vector<unsigned> indices_;
  std::vector<unsigned> section_;
  std::map<DetectorEnum,unsigned> enumMap_;

  unsigned nLayers_;
  unsigned nSections_;
  bool bypassRadius_;
};

HGCSSDetector & theDetector();




#endif
