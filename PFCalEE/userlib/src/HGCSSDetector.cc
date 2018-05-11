#include "HGCSSDetector.hh"
#include <iostream>

HGCSSDetector & theDetector(){
  static HGCSSDetector lDet;
  static bool firstDet=true;
  if (firstDet) std::cout << " -- Created detector static object." << std::endl;
  firstDet=false;
  return lDet;
}

//
void HGCSSDetector::buildDetector( HGCSSInfo * info) {

  using namespace std;

  bypassRadius_ = false;  
  reset();

  //
  size_t nsens(info->sensitiveZ().size());
  indices_.clear();
  indices_.resize(7,0);
  indices_[0] = 0;
  indices_[1] = 10;
  indices_[2] = 20;
  indices_[3] = 28;
  indices_[4] = 53;
  indices_[5] = 57;
  indices_[6] = nsens;
  sensitiveZ_.resize(nsens,0);
  etaBoundary_.resize(nsens,0);

  int transLayer(52);
  if(info->version()==65 || info->version()==650 || info->version()==66 || info->version()==660) {
    transLayer  = 46;
    indices_[4] = 47;
    indices_[5] = 51;
  }

  for(size_t i=0; i<nsens; i++) {
    sensitiveZ_[i]  = info->sensitiveZ()[i];
    etaBoundary_[i] = info->etaBoundaryMin()[i];    
    if(info->version()>600 && i<transLayer) { etaBoundary_[i]=0; sensitiveZ_[i]=0; }
    if(i==transLayer) etaBoundary_[i]=0;
    if(i>transLayer)  etaBoundary_[i]=info->etaBoundaryMax()[i];
  }
  
  initSubDetectors(false,info->version());
  finishInitialisation();
}

//
void HGCSSDetector::buildDetector(const unsigned versionNumber,
				  bool concept,
				  bool isCaliceHcal,
				  bool bypassR){
  
  bypassRadius_ = bypassR;
  reset();
  initialiseIndices(versionNumber);
  initSubDetectors(isCaliceHcal,versionNumber);
  finishInitialisation();
}

//
void HGCSSDetector::initSubDetectors(bool isCaliceHcal,int versionNumber) { 
  HGCSSSubDetector FECAL;
  FECAL.type = DetectorEnum::FECAL;
  FECAL.name = "FECAL";
  FECAL.layerIdMin = indices_[0];
  FECAL.layerIdMax = indices_[1];
  FECAL.mipWeight = 1./0.0548;//mip for 200um si
  FECAL.absWeight = 1.;//ratio of abs dedx
  FECAL.gevWeight = 1.0;
  FECAL.gevOffset = 0.0;
  FECAL.isSi = true;
  //if (versionNumber>=30) 
  if (!bypassRadius_) FECAL.radiusLim = 750;
  else FECAL.radiusLim = 0;
  if (FECAL.nLayers()>0) theDetector().addSubdetector(FECAL);
  
  HGCSSSubDetector MECAL;
  MECAL.type = DetectorEnum::MECAL;
  MECAL.name = "MECAL";
  MECAL.layerIdMin = indices_[1];
  MECAL.layerIdMax = indices_[2];
  MECAL.mipWeight = 1./0.0548;//mip for 200um si
  MECAL.absWeight = 8.001/5.848;//ratio of abs dedx
  MECAL.gevWeight = 1.0;
  MECAL.gevOffset = 0.0;
  MECAL.isSi = true;
  //if (versionNumber>=30) 
  if (!bypassRadius_) MECAL.radiusLim = 750;
  else MECAL.radiusLim = 0;
  if (MECAL.nLayers()>0) theDetector().addSubdetector(MECAL);
  
  HGCSSSubDetector BECAL;
  BECAL.type = DetectorEnum::BECAL;
  BECAL.name = "BECAL";
  BECAL.layerIdMin = indices_[2];
  BECAL.layerIdMax = indices_[3];
  BECAL.mipWeight = 1./0.0548;//mip for 200um si
  BECAL.absWeight = 10.854/5.848;//ratio of abs dedx
  BECAL.gevWeight = 1.0;
  BECAL.gevOffset = 0.0;
  BECAL.isSi = true;
  //if (versionNumber>=30) 
  if (!bypassRadius_) BECAL.radiusLim = 750;
  else BECAL.radiusLim = 0;
  if (BECAL.nLayers()>0) theDetector().addSubdetector(BECAL);
  
  HGCSSSubDetector FHCAL;
  FHCAL.type = DetectorEnum::FHCAL;
  FHCAL.name = "FHCAL";
  FHCAL.layerIdMin = indices_[3];
  FHCAL.layerIdMax = indices_[4];
  FHCAL.mipWeight = 1./0.0548;//mip for 200um si
  FHCAL.absWeight = 65.235/5.848;//ratio of abs dedx
  //if (!concept) FHCAL.absWeight = 0.5*65.235/5.848;
  FHCAL.gevWeight = 1.;
  FHCAL.gevOffset = 0.;
  FHCAL.isSi = true;
  //if (versionNumber>=30) 
  if (!bypassRadius_) FHCAL.radiusLim = 600;
  else FHCAL.radiusLim = 0;
  if (isCaliceHcal) {
    FHCAL.mipWeight = 1./0.816;
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
  BHCAL.layerIdMin = indices_[4];
  BHCAL.layerIdMax = indices_[5];
  BHCAL.mipWeight = 1./0.63;//was 1.49 for 9mm scint
  if ((versionNumber>59 && versionNumber<67) || (versionNumber>=630 && versionNumber<=660)) 
    BHCAL.mipWeight = 1./0.497;//for 3mm scint, v8
  BHCAL.absWeight = 1.0;//92.196/5.848;
  BHCAL.gevWeight = 1.0;
  BHCAL.gevOffset = 0.0;
  BHCAL.isScint = true;
  
  if (isCaliceHcal) {
    BHCAL.name = "BHCAL1";
    BHCAL.mipWeight = FHCAL.mipWeight;
    BHCAL.absWeight = 1.;
    BHCAL.gevWeight = FHCAL.gevWeight;//MIPtoGeV
    BHCAL.gevOffset = 0.0;
  }
  if (BHCAL.nLayers()>0) theDetector().addSubdetector(BHCAL);
  
  HGCSSSubDetector BHCAL2;
  BHCAL2.type = DetectorEnum::BHCAL2;
  BHCAL2.name = "BHCAL2";
  BHCAL2.layerIdMin = indices_[5];
  BHCAL2.layerIdMax = indices_[6];
  BHCAL2.mipWeight = BHCAL.mipWeight;
  BHCAL2.absWeight = 1;
  BHCAL2.gevWeight = 1;

  if (isCaliceHcal) {
    BHCAL2.mipWeight = FHCAL.mipWeight;
    BHCAL2.absWeight = 104./21.;
    BHCAL2.gevWeight = FHCAL.gevWeight;//MIPtoGeV
  }
  BHCAL2.gevOffset = 0.0;
  BHCAL2.isScint = true;
  
  if (BHCAL2.nLayers()>0) theDetector().addSubdetector(BHCAL2);
}

//
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
  //indices_.push_back(adet.layerIdMin);
}
  
void HGCSSDetector::finishInitialisation(){
  nSections_ = subdets_.size();
  //indices_.push_back(subdets_[nSections_-1].layerIdMax);
  const unsigned lastEle = indices_.size()-1;
  nLayers_ = indices_[lastEle];

  unsigned secIndex[lastEle];
  unsigned iS(0);
  for(unsigned i(0); i<lastEle ;i++){
     secIndex[i] = iS;
     if(indices_[i] < indices_[i+1])iS +=1;
  }
  //initialise layer-section conversion
  section_.resize(nLayers_,0);
  for (unsigned iL(0); iL<nLayers_;++iL){
    for (unsigned i(0); i<lastEle;++i){
      if (iL >= indices_[i] && iL < indices_[i+1]) section_[iL] = secIndex[i];
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
