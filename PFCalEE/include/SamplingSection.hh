#ifndef _samplingsection_hh_
#define _samplingsection_hh_

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Colour.hh"

#include <iomanip>
#include <vector>

#include "G4SiHit.hh"

class SamplingSection
{
public:
  //CTOR
  SamplingSection(const std::vector<G4double> & aThicknessVec,
		  const std::vector<std::string> & aMaterialVec)
    : isCylinder_(false),
      radius_(-1),
      fullySensitive_(false)
  {
    if (aMaterialVec.size() != aThicknessVec.size()) {
      G4cout << " -- ERROR in sampling section definition. Expect input vectors with same size containing thicknesses (size=" << aThicknessVec.size() << ") and material names (size=" << aMaterialVec.size() << "). Exiting..." << G4endl;
      exit(1);
    }
    Total_thick = 0;
    n_sens_elements=0;
    n_elements=0;
    n_sectors=0;
    ele_thick.clear();
    ele_name.clear();
    ele_X0.clear();
    ele_L0.clear();
    ele_vol.clear();
    for (unsigned ie(0);  ie<aThicknessVec.size(); ++ie){
      //consider only material with some non-0 width...
      if (aThicknessVec[ie]>0){
	ele_thick.push_back(aThicknessVec[ie]);
	ele_name.push_back(aMaterialVec[ie]);
	ele_X0.push_back(0);
	ele_dEdx.push_back(0);
	ele_L0.push_back(0);
	ele_vol.push_back(0);
	Total_thick+=aThicknessVec[ie];
	++n_elements;
	//the following method check the total size...
	//so incrementing first.
	if (isSensitiveElement(n_elements-1)) {
	  G4SiHitVec lVec;
	  sens_HitVec.push_back(lVec);
	  ++n_sens_elements;
	}
      }
    }
    sens_HitVec_size_max = 0;
    resetCounters();

    std::cout << " -- End of sampling section initialisation. Input " << aThicknessVec.size() << " elements, constructing " << n_elements << " elements with " << n_sens_elements << " sensitive elements." << std::endl;

  };

  bool isCylinder_;
  float radius_;
  void setSectionAsCylinder(float radius) { isCylinder_=true; radius_=radius; }
  bool fullySensitive_;
  void enableAsSensitive() { fullySensitive_=true;}

  //DTOR
  ~SamplingSection() { };

  inline void setNumberOfSectors(const unsigned nSec){
    n_sectors = nSec;
    ele_vol.clear();
    for (unsigned ie(0);  ie<n_elements*n_sectors; ++ie){
      ele_vol.push_back(0);
    }
  };

  //
  void add(G4double den, G4double dl, 
	   G4double globalTime,G4int pdgId,
	   G4VPhysicalVolume* vol, 
	   const G4ThreeVector & position,
	   G4int trackID, G4int parentID,
	   G4int layerId);
  
  inline bool isSensitiveElement(const unsigned & aEle){
    if (aEle < n_elements &&
	(ele_name[aEle] == "Si" || ele_name[aEle] == "Scintillator")
	) return true;
    return false;
  };

  inline unsigned getSensitiveLayerIndex(std::string astr){
    if (astr.find("_")== astr.npos) return 0;
    size_t pos = astr.find("phys");
    //std::cout << astr << " " << pos << std::endl;
    if (pos != astr.npos && pos>1) {
      unsigned idx = 0;//atoi(astr[pos-1]);
      std::istringstream(astr.substr(pos-1,1))>>idx;
      return idx;
    }
    return 0;
  };

  inline G4Colour g4Colour(const unsigned & aEle){
    if (isSensitiveElement(aEle)) return G4Colour::Red();
    if (ele_name[aEle] == "Cu") return G4Colour::Black();
    if (isAbsorberElement(aEle)) return G4Colour::Gray();
    if (ele_name[aEle] == "PCB") return G4Colour::Blue();
    if (ele_name[aEle] == "Air") return G4Colour::Cyan();
    return G4Colour::Yellow();
  };

  inline bool isAbsorberElement(const unsigned & aEle){
    if (aEle < n_elements && 
	(
	 ele_name[aEle] == "Pb" || ele_name[aEle] == "Cu" || 
	 ele_name[aEle] == "W" || ele_name[aEle] == "Brass" ||
	 ele_name[aEle] == "Fe" || ele_name[aEle] == "Steel" || 
	 ele_name[aEle] == "SSteel" || ele_name[aEle] == "Al" ||
	 ele_name[aEle] == "WCu" || ele_name[aEle] == "NeutMod"
	 )
	) return true;
    return false;
  };

  //reset
  inline void resetCounters()
  {
    ele_den.resize(n_elements,0);
    ele_dl.resize(n_elements,0);
    sens_time.resize(n_sens_elements,0);
    sens_gFlux.resize(n_sens_elements,0);
    sens_eFlux.resize(n_sens_elements,0);
    sens_muFlux.resize(n_sens_elements,0);
    sens_neutronFlux.resize(n_sens_elements,0);
    sens_hadFlux.resize(n_sens_elements,0);
    //reserve some space based on first event....
    for (unsigned idx(0); idx<n_sens_elements; ++idx){
      if (sens_HitVec[idx].size() > sens_HitVec_size_max) {
	sens_HitVec_size_max = 2*sens_HitVec[idx].size();
	G4cout << "-- SamplingSection::resetCounters(), space reserved for HitVec vector increased to " << sens_HitVec_size_max << G4endl;
      }
      sens_HitVec[idx].clear();
      sens_HitVec[idx].reserve(sens_HitVec_size_max);
    }
  }
  
  //
  G4double getMeasuredEnergy(bool weighted=true);
  G4double getAbsorbedEnergy();
  G4double getTotalEnergy();
  G4double getAbsorberX0();  
  G4double getAbsorberdEdx();  
  G4double getAbsorberLambda();
  G4double getHadronicFraction();
  G4double getNeutronFraction();
  G4double getMuonFraction();
  G4double getPhotonFraction();
  G4double getElectronFraction();
  G4double getAverageTime();
  G4int getTotalSensHits();
  G4double getTotalSensE();

  const G4SiHitVec & getSiHitVec(const unsigned & idx) const;
  void trackParticleHistory(const unsigned & idx,const G4SiHitVec & incoming);

  //
  void report(bool header=false);

  //members
  unsigned n_elements;
  unsigned n_sectors;
  unsigned n_sens_elements;
  std::vector<G4double>           ele_thick;
  std::vector<std::string>        ele_name;
  std::vector<G4double>           ele_X0;
  std::vector<G4double>           ele_dEdx;
  std::vector<G4double>           ele_L0;
  std::vector<G4double>           ele_den;
  std::vector<G4double>           ele_dl;
  std::vector<G4VPhysicalVolume*> ele_vol;
  std::vector<G4double>           sens_gFlux, sens_eFlux, sens_muFlux, sens_neutronFlux, sens_hadFlux, sens_time;
  G4double Total_thick;
  std::vector<G4SiHitVec> sens_HitVec;
  unsigned sens_HitVec_size_max;


};

#endif
