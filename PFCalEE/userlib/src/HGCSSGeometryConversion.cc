#include "HGCSSGeometryConversion.hh"
#include <sstream>
#include <iostream>
#include <cmath>

HGCSSGeometryConversion::HGCSSGeometryConversion(std::string filePath, int model, double cellsize){

  width_ = 200;//mm
  model_ = model;
  if (model==1) width_ = 500;
  else if (model == 2) width_ = 1700*2;
  else if (model == 3) width_ = 1000;
  
  cellSize_ = cellsize;
  
}

unsigned HGCSSGeometryConversion::getNumberOfSiLayers(const DetectorEnum type,
						      const double & eta) const
{
  if (model_ != 2) return 3;
  if (type == DetectorEnum::FHCAL) return 3;
  unsigned etaBin = 0;
  if (fabs(eta)>=1.4 && fabs(eta)<1.9) etaBin = 1;
  else if (fabs(eta)>=1.9 && fabs(eta)<2.4) etaBin = 2;
  else if (fabs(eta) >= 2.4) etaBin = 3;
  if (etaBin==0){
    if (type == DetectorEnum::FECAL) return 2;
    else if (type == DetectorEnum::MECAL) return 2;
    else if (type == DetectorEnum::BECAL) return 2;
  }
  else {
    if (etaBin==1) return 3;
    else if (etaBin==2) return 2;
    else return 1;
  }
  return 3;
}

HGCSSGeometryConversion::~HGCSSGeometryConversion(){
  std::map<DetectorEnum,std::vector<TH2D *> >::iterator liter =
    HistMapE_.begin();
  for (; liter !=HistMapE_.end();++liter){
    deleteHistos(liter->second);
  }
  HistMapE_.clear();
  liter = HistMapTime_.begin();
  for (; liter !=HistMapTime_.end();++liter){
    deleteHistos(liter->second);
  }
  HistMapTime_.clear();
  liter = HistMapZ_.begin();
  for (; liter !=HistMapZ_.end();++liter){
    deleteHistos(liter->second);
  }
  HistMapZ_.clear();
}


void HGCSSGeometryConversion::deleteHistos(std::vector<TH2D *> & aVec){
  if (aVec.size()!=0){
    for (unsigned iL(0); iL<aVec.size();++iL){
      aVec[iL]->Delete();
    }
  }
}

void HGCSSGeometryConversion::setGranularity(const std::vector<unsigned> & granul){
  granularity_.reserve(granul.size());
  for (unsigned iL(0); iL<granul.size();++iL){
    granularity_.push_back(granul[iL]);
  }
}


 void HGCSSGeometryConversion::initialiseHistos(const bool recreate){

   for (unsigned iS(0); iS<theDetector().nSections();++iS){
     resetVector(HistMapE_[theDetector().detType(iS)],"EmipHits",theDetector().detName(iS),theDetector().subDetectorBySection(iS),theDetector().nLayers(iS),recreate);
     //std::cout << " check: " << HistMapE_[theDetector().detType(iS)].size() << std::endl;
     
     std::vector<double> avgvecE;
     avgvecE.resize(theDetector().nLayers(iS),0);
     avgMapE_[theDetector().detType(iS)]=avgvecE;
     
     resetVector(HistMapTime_[theDetector().detType(iS)],"TimeHits",theDetector().detName(iS),theDetector().subDetectorBySection(iS),theDetector().nLayers(iS),recreate);
     //std::cout << " check: " << HistMapTime_[theDetector().detType(iS)].size() << std::endl;

     resetVector(HistMapZ_[theDetector().detType(iS)],"zHits",theDetector().detName(iS),theDetector().subDetectorBySection(iS),theDetector().nLayers(iS),recreate);
     //std::cout << " check: " << HistMapZ_[theDetector().detType(iS)].size() << std::endl;
     
     std::vector<double> avgvecZ;
     avgvecZ.resize(theDetector().nLayers(iS),0);
     avgMapZ_[theDetector().detType(iS)]=avgvecZ;
   }
 }

void HGCSSGeometryConversion::fill(const DetectorEnum type,
				   const unsigned newlayer,
				   const double & weightedE,
				   const double & aTime,
				   const double & posx,
				   const double & posy,
				   const double & posz)
{
  HistMapE_[type][newlayer]->Fill(posx,posy,weightedE);
  HistMapTime_[type][newlayer]->Fill(posx,posy,weightedE*aTime);
  HistMapZ_[type][newlayer]->Fill(posx,posy,weightedE*posz);
  avgMapZ_[type][newlayer] += weightedE*posz;
  avgMapE_[type][newlayer] += weightedE;
}

double HGCSSGeometryConversion::getAverageZ(const unsigned layer){
  const HGCSSSubDetector & subdet = theDetector().subDetectorByLayer(layer);
  unsigned newlayer = layer-subdet.layerIdMin;
  double avg = 0;
  if (avgMapE_[subdet.type][newlayer]>0)
    avg =avgMapZ_[subdet.type][newlayer]/avgMapE_[subdet.type][newlayer];
  return avg;
}

TH2D * HGCSSGeometryConversion::get2DHist(const unsigned layer,std::string name){
  const HGCSSSubDetector & subdet = theDetector().subDetectorByLayer(layer);
  unsigned newlayer = layer-subdet.layerIdMin;
  if (name == "E") return HistMapE_[subdet.type][newlayer];
  else if (name == "Time") return HistMapTime_[subdet.type][newlayer];
  else if (name == "Z") return HistMapZ_[subdet.type][newlayer];
  else {
    std::cerr << " ERROR !! Unknown histogram name. Exiting..." << std::endl;
    exit(1);
  }
}

unsigned HGCSSGeometryConversion::getGranularity(const unsigned aLayer, const HGCSSSubDetector & adet){
  unsigned idx = adet.layerIdMin+aLayer;
  return granularity_[idx];
}

void HGCSSGeometryConversion::resetVector(std::vector<TH2D *> & aVec,
					  std::string aVar,
					  std::string aString,
					  const HGCSSSubDetector & aDet,
					  const unsigned nLayers,
					  bool recreate)
{
  //std::cout << " vector size: " << aVar << " " << aString << " = " << aVec.size() << std::endl;
  if (recreate){
    for (unsigned iL(0); iL<aVec.size();++iL){
      aVec[iL]->Delete();
    }
    aVec.clear();
  }
  if (aVec.size()!=0){
    for (unsigned iL(0); iL<aVec.size();++iL){
      aVec[iL]->Reset();
    }
  }
  else {
    if (nLayers > 0){
      aVec.resize(nLayers,0);
      std::cout << " -- Creating " << nLayers << " 2D histograms for " << aVar << " " << aString 
	//<< " with " << nBins << " bins between " << min << " and " << max 
		<< std::endl;
      for (unsigned iL(0); iL<nLayers;++iL){
	std::ostringstream lname;
	lname << aVar << "_" << aString << "_" << iL ;
	double newcellsize = cellSize_*getGranularity(iL,aDet);
	//take smallest pair integer to be sure to fit in the geometry
	//even if small dead area introduced at the edge
	unsigned nBins = static_cast<unsigned>(width_*1./(newcellsize*2.))*2;
	double min = -1.0*nBins*newcellsize/2.;
	double max = nBins*newcellsize/2.;
	aVec[iL] = new TH2D(lname.str().c_str(),";x(mm);y(mm)",nBins,min,max,nBins,min,max);
	if (iL==0) std::cout << " ---- bins, min, max = " << nBins << " " << min << " " << max << std::endl;
      }
    }
  }
  //std::cout << " vector size after: " << aVar << " " << aString << " = " << aVec.size() << std::endl;
}





