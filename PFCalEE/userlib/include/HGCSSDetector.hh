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
    isScint(false)
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

  inline unsigned nLayers() const{
    return (layerIdMax-layerIdMin);
  };

private:

};

class HGCSSDetector {

public:
  friend HGCSSDetector & theDetector();

  void buildDetector(const unsigned indices[],
		     bool concept=true,
		     bool isCaliceHcal=false);

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
  HGCSSDetector(){};

  ~HGCSSDetector(){
    reset();
  };
  
  std::vector<HGCSSSubDetector> subdets_;
  std::vector<unsigned> indices_;
  std::vector<unsigned> section_;
  std::map<DetectorEnum,unsigned> enumMap_;

  unsigned nLayers_;
  unsigned nSections_;
};

HGCSSDetector & theDetector();




#endif
