#include "G4VPhysicalVolume.hh"

#include "SamplingSection.hh"

//
void SamplingSection::add(G4double den, G4double dl, 
			  G4double globalTime, G4int pdgId, 
			  G4VPhysicalVolume* vol, 
			  const G4ThreeVector & position,
			  G4int trackID, G4int parentID,
			  G4int layerId)
{
  std::string lstr = vol->GetName();
  for (unsigned ie(0); ie<n_elements;++ie){
    if(ele_vol[ie] && lstr==ele_vol[ie]->GetName()){ 
      ele_den[ie]+=den;
      ele_dl[ie]+=dl;
      if (isSensitiveElement(ie)) {//if Si || sci
	unsigned idx = getSensitiveLayerIndex(lstr);
	//std::cout << "sens layer " << idx << " " << lstr << std::endl;
	sens_time[idx]+=den*globalTime;
	
	//discriminate further by particle type
	if(abs(pdgId)==22)      sens_gFlux[idx] += den;
	else if(abs(pdgId)==11) sens_eFlux[idx] += den;
	else if(abs(pdgId)==13) sens_muFlux[idx] += den;
	else if (abs(pdgId)==2112) sens_neutronFlux[idx] += den;
	else {
	  sens_hadFlux[idx] += den;
	}
	
	//add hit
	G4SiHit lHit;
	lHit.energy = den;
	lHit.time = globalTime;
	lHit.pdgId = pdgId;
	lHit.layer = layerId;
	lHit.hit_x = position.x();
	lHit.hit_y = position.y();
	lHit.hit_z = position.z();
	lHit.trackId = trackID;
	lHit.parentId = parentID;
	sens_HitVec[idx].push_back(lHit);

      }//if Si
    }//if in right material
    
  }//loop on available materials
  
}

//
void SamplingSection::report(bool header)
{
  if(header) G4cout << "E/[MeV]\t  Si\tAbsorber\tTotal\tSi g frac\tSi e frac\tSi mu frac\tSi had frac\tSi <t> \t nG4SiHits" << G4endl;
  G4cout << std::setprecision(3) << "\t  " << getMeasuredEnergy(false) << "\t" << getAbsorbedEnergy() << "\t\t" << getTotalEnergy() << "\t"
	 << getPhotonFraction() << "\t" << getElectronFraction() << "\t" << getMuonFraction() << "\t" << getHadronicFraction() << "\t"
	 << getAverageTime() << "\t"
	 << getTotalSensHits() << "\t"
	 << G4endl; 
}

G4int SamplingSection::getTotalSensHits()
{
  G4int tot=0;
  for (unsigned ie(0); ie<n_sens_elements;++ie){
    tot += sens_HitVec[ie].size();
  }
  return tot;
}

G4double SamplingSection::getTotalSensE(){
  double etot = 0;
  for (unsigned ie(0); ie<n_elements;++ie) {
    if (isSensitiveElement(ie)) etot += ele_den[ie];
  }
  return etot;
}

G4double SamplingSection::getAverageTime()
{
  double etot = getTotalSensE();
  double time = 0;
  for (unsigned ie(0); ie<n_sens_elements;++ie){
    time += sens_time[ie];
  }
  return etot>0 ? time/etot : 0 ;
}

//
G4double SamplingSection::getPhotonFraction()
{
  double etot = getTotalSensE();
  double val = 0;
  for (unsigned ie(0); ie<n_sens_elements;++ie){
    val += sens_gFlux[ie];
  }
  return etot>0 ? val/etot : 0 ;
}

//
G4double SamplingSection::getElectronFraction()
{
  double etot = getTotalSensE();
  double val = 0;
  for (unsigned ie(0); ie<n_sens_elements;++ie){
    val += sens_eFlux[ie];
  }
  return etot>0 ? val/etot : 0 ;
}

//
G4double SamplingSection::getMuonFraction()
{
  double etot = getTotalSensE();
  double val = 0;
  for (unsigned ie(0); ie<n_sens_elements;++ie){
    val += sens_muFlux[ie];
  }
  return etot>0 ? val/etot : 0 ;
}

//
G4double SamplingSection::getNeutronFraction()
{
  double etot = getTotalSensE();
  double val = 0;
  for (unsigned ie(0); ie<n_sens_elements;++ie){
    val += sens_neutronFlux[ie];
  }
  return etot>0 ? val/etot : 0 ;
}

//
G4double SamplingSection::getHadronicFraction()
{
  double etot = getTotalSensE();
  double val = 0;
  for (unsigned ie(0); ie<n_sens_elements;++ie){
    val += sens_hadFlux[ie];
  }
  return etot>0 ? val/etot : 0 ;
}

//
G4double SamplingSection::getMeasuredEnergy(bool weighted)
{
  G4double weight=(weighted ? getAbsorberX0() : 1.0);
  return weight*getTotalSensE();
}

//
G4double SamplingSection::getAbsorberX0()
{
  // G4cout << Pb_thick << " " << Pb_X0 << G4endl;
  double val=0;  
  for (unsigned ie(0); ie<n_elements;++ie){
    if (isAbsorberElement(ie))
      if (ele_X0[ie]>0) val += ele_thick[ie]/ele_X0[ie];
  }
  return val;
}

//
G4double SamplingSection::getAbsorberLambda()
{
  double val=0;  
  for (unsigned ie(0); ie<n_elements;++ie){
    if (isAbsorberElement(ie) && ele_L0[ie]>0)
      val += ele_thick[ie]/ele_L0[ie];
  }
  return val;
}

//
G4double SamplingSection::getAbsorbedEnergy()
{
  double val=0;  
  for (unsigned ie(0); ie<n_elements;++ie){
    if (isAbsorberElement(ie)) val += ele_den[ie];
  }
  return val;
}

//
G4double SamplingSection::getTotalEnergy()
{
 double val=0;  
  for (unsigned ie(0); ie<n_elements;++ie){
    val += ele_den[ie];
  }
  return val;
}

const G4SiHitVec & SamplingSection::getSiHitVec(const unsigned & idx) const
{
  return sens_HitVec[idx];
}

void SamplingSection::trackParticleHistory(const unsigned & idx, const G4SiHitVec & incoming)
{
  for (unsigned iP(0); iP<sens_HitVec[idx].size(); ++iP){//loop on g4hits
    G4int parId = sens_HitVec[idx][iP].parentId;
    for (unsigned iI(0); iI<incoming.size(); ++iI){//loop on previous layer
      G4int trId = incoming[iI].trackId;
      if (trId == parId) sens_HitVec[idx][iP].parentId = incoming[iI].parentId;
    }//loop on previous layer
  }//loop on g4hits
  //  return sens_HitVec;
}
