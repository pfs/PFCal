#include "HGCSSDetector.hh"
#include <iostream>

HGCSSDetector & theDetector(){
  static HGCSSDetector lDet;
  static bool firstDet=true;
  if (firstDet) std::cout << " -- Created detector static object." << std::endl;
  firstDet=false;
  return lDet;
}

void HGCSSDetector::buildDetector(const unsigned indices[],
				  bool concept,
				  bool isCaliceHcal){
  
  reset();
  HGCSSSubDetector FECAL;
  FECAL.type = DetectorEnum::FECAL;
  FECAL.name = "FECAL";
  FECAL.layerIdMin = indices[0];
  FECAL.layerIdMax = indices[1];
  FECAL.mipWeight = 2./0.0548;//mip for 100um si
  FECAL.absWeight = 1.;//ratio of abs dedx
  FECAL.gevWeight = 1.0;
  FECAL.gevOffset = 0.0;
  FECAL.isSi = true;
  if (FECAL.nLayers()>0) theDetector().addSubdetector(FECAL);
  
  HGCSSSubDetector MECAL;
  MECAL.type = DetectorEnum::MECAL;
  MECAL.name = "MECAL";
  MECAL.layerIdMin = indices[1];
  MECAL.layerIdMax = indices[2];
  MECAL.mipWeight = 2./0.0548;//mip for 100um si
  MECAL.absWeight = 8.001/5.848;//ratio of abs dedx
  MECAL.gevWeight = 1.0;
  MECAL.gevOffset = 0.0;
  MECAL.isSi = true;
  if (MECAL.nLayers()>0) theDetector().addSubdetector(MECAL);
  
  HGCSSSubDetector BECAL;
  BECAL.type = DetectorEnum::BECAL;
  BECAL.name = "BECAL";
  BECAL.layerIdMin = indices[2];
  BECAL.layerIdMax = indices[3];
  BECAL.mipWeight = 2./0.0548;//mip for 100um si
  BECAL.absWeight = 10.854/5.848;//ratio of abs dedx
  BECAL.gevWeight = 1.0;
  BECAL.gevOffset = 0.0;
  BECAL.isSi = true;
  if (BECAL.nLayers()>0) theDetector().addSubdetector(BECAL);
  
  HGCSSSubDetector FHCAL;
  FHCAL.type = DetectorEnum::FHCAL;
  FHCAL.name = "FHCAL";
  FHCAL.layerIdMin = indices[3];
  FHCAL.layerIdMax = indices[4];
  FHCAL.mipWeight = 3./0.0849;//mip for 100um si
  FHCAL.absWeight = 65.235/5.848;//ratio of abs dedx
  if (!concept) FHCAL.absWeight = 0.5*65.235/5.848;
  FHCAL.gevWeight = 1.;
  FHCAL.gevOffset = 0.;
  FHCAL.isSi = true;
  if (isCaliceHcal) {
    FHCAL.mipWeight = 1./0.807;
    FHCAL.absWeight = 1.;
    FHCAL.gevWeight = 1./39.32;//MIPtoGeV
    FHCAL.gevOffset = -1.8/39.32;//offset in GeV
    FHCAL.isScint = true;
    FHCAL.isSi = false;
  }
  if (FHCAL.nLayers()>0) theDetector().addSubdetector(FHCAL);
  
  HGCSSSubDetector BHCAL;
  BHCAL.type = DetectorEnum::BHCAL1;
  BHCAL.name = "BHCAL";
  BHCAL.layerIdMin = indices[4];
  BHCAL.layerIdMax = indices[5];
  BHCAL.mipWeight = 1./1.49;
  BHCAL.absWeight = 92.196/5.848;
  BHCAL.gevWeight = 1.0;
  BHCAL.gevOffset = 0.0;
  BHCAL.isScint = true;
  
  if (isCaliceHcal) {
    BHCAL.name = "BHCAL1";
    BHCAL.mipWeight = 1./0.807;
    BHCAL.absWeight = 1.;
    BHCAL.gevWeight = FHCAL.gevWeight;//MIPtoGeV
    BHCAL.gevOffset = 0.0;
  }
  if (BHCAL.nLayers()>0) theDetector().addSubdetector(BHCAL);
  
  HGCSSSubDetector BHCAL2;
  BHCAL2.type = DetectorEnum::BHCAL2;
  BHCAL2.name = "BHCAL2";
  BHCAL2.layerIdMin = indices[5];
  BHCAL2.layerIdMax = indices[6];
  BHCAL2.mipWeight = 1./0.807;
  BHCAL2.absWeight = 104./21.;
  BHCAL2.gevWeight = FHCAL.gevWeight;//MIPtoGeV
  BHCAL2.gevOffset = 0.0;
  BHCAL2.isScint = true;
  
  if (BHCAL2.nLayers()>0) theDetector().addSubdetector(BHCAL2);
  
  finishInitialisation();
  
}


const HGCSSSubDetector & HGCSSDetector::subDetectorByLayer(const unsigned aLayer){
  unsigned section = getSection(aLayer);
  return subdets_[section];
}

unsigned HGCSSDetector::getSection(const unsigned aLayer) const{
  if (aLayer>=nLayers_) {
    std::cerr << " -- Error ! Trying to access layer " << aLayer 
	      << " outside of range. nLayers = " << nLayers_
	      << std::endl;
    exit(1);
  }
  return section_[aLayer];
}

void HGCSSDetector::addSubdetector(const HGCSSSubDetector & adet){
  subdets_.push_back(adet);
  enumMap_[adet.type]=subdets_.size()-1;
  indices_.push_back(adet.layerIdMin);
}
  
void HGCSSDetector::finishInitialisation(){
  nSections_ = subdets_.size();
  indices_.push_back(subdets_[nSections_-1].layerIdMax);
  unsigned lastEle = indices_.size()-1;
  nLayers_ = indices_[lastEle];
  //initialise layer-section conversion
  section_.resize(nLayers_,0);
  for (unsigned iL(0); iL<nLayers_;++iL){
    for (unsigned i(0); i<lastEle;++i){
      if (iL >= indices_[i] && iL < indices_[i+1]) section_[iL] = i;
    }
  }
  printDetector(std::cout);
}

const HGCSSSubDetector & HGCSSDetector::subDetectorByEnum(DetectorEnum adet){
  if (enumMap_.find(adet) == enumMap_.end()){
    std::cerr << " -- Error ! Trying to access subdetector enum not present in this detector: "
	      << adet 
	      << std::endl;
    exit(1);
  } 
  return subdets_[enumMap_[adet]];
}

void HGCSSDetector::reset() {
  subdets_.clear();
  enumMap_.clear();
  indices_.clear();
  section_.clear();
}

void HGCSSDetector::printDetector(std::ostream & aOs) const{
  std::cout << " -------------------------- " << std::endl
	    << " -- Detector information -- " << std::endl
	    << " -------------------------- " << std::endl
	    << " - nSections = " << nSections_ << std::endl
	    << " - nLayers = " << nLayers_ << std::endl
	    << " - detNames = " ;
  for (unsigned i(0); i<nSections_;++i){
    std::cout << " " << detName(i);
  }
  std::cout << std::endl;
  std::cout << " -------------------------- " << std::endl;
}
