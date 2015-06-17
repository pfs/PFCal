#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include <boost/algorithm/string.hpp>

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4CSGSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4PhysicalConstants.hh"

using namespace std;

//
DetectorConstruction::DetectorConstruction(G4int ver, G4int mod,
					   std::string absThickW,
					   std::string absThickPb,
					   std::string dropLayer) : 
  version_(ver), model_(mod), addPrePCB_(false)
{
  SetWThick(absThickW);
  SetPbThick(absThickPb);
  SetDropLayers(dropLayer);

  //radiation lengths: cf. http://pdg.lbl.gov/2012/AtomicNuclearProperties/
  //W 3.504 mm
  //Pb 5.612 mm
  //Cu 14.36 mm
  switch(version_)
    {
      //cf. http://arxiv.org/abs/0805.4833
    case v_CALICE:
      {
	G4cout << "[DetectorConstruction] starting v_CALICE (10x0.4+10x0.8+10x1.2)X_0 with Tungsten" << G4endl;
	std::vector<G4double> lThick;
	std::vector<std::string> lEle;
	lThick.push_back(0.4*3.504*mm);lEle.push_back("W");
	lThick.push_back(0.525*mm);lEle.push_back("Si");
	lThick.push_back(1.0*mm);lEle.push_back("PCB");
	lThick.push_back(2.5*mm);lEle.push_back("Air");

	for(unsigned i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	lThick[0] = 0.8*3.504*mm;
	for(unsigned i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	lThick[0] = 1.2*3.504*mm;
	for(unsigned i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	break;
      }

    case v_HGCALEE_v6: case v_HGCAL_v6:
      {
	G4cout << "[DetectorConstruction] starting v_HCALEE_v6"<< G4endl;
	G4double airThick = 3*mm;
	G4double pcbThick = 2*mm;

  std::vector<G4double> lThickL;
  std::vector<std::string> lEleL;
  // Moderator
  lThickL.push_back(180*mm);lEleL.push_back("NeutMod");
  // Thermal insulator
  lThickL.push_back(2.*mm);lEleL.push_back("Al");
  lThickL.push_back(26.*mm);lEleL.push_back("Foam");
  lThickL.push_back(2.*mm);lEleL.push_back("Al");
  // Pre-shower
  lThickL.push_back(0.1*mm);lEleL.push_back("Si");
  lThickL.push_back(0.1*mm);lEleL.push_back("Si");
  lThickL.push_back(0.1*mm);lEleL.push_back("Si");
  m_caloStruct.push_back( SamplingSection(lThickL,lEleL) );

  lThickL.clear();
  lEleL.clear();
  unsigned nFEModsNom = 5;
  unsigned nFEModsSub = 0;
  unsigned nFEModsNet = nFEModsNom - nFEModsSub;
  double rScale = ((double) nFEModsNom )/((double) nFEModsNet);
  // Front module
  lThickL.push_back(0.5*mm*rScale);lEleL.push_back("Cu");
  lThickL.push_back(0.5*mm);lEleL.push_back("CFMix");
  lThickL.push_back(2.0*mm*rScale);lEleL.push_back("W");
  lThickL.push_back(0.5*mm);lEleL.push_back("CFMix");
  lThickL.push_back(0.5*mm*rScale);lEleL.push_back("Cu");
  lThickL.push_back(airThick);lEleL.push_back("Air");
  lThickL.push_back(pcbThick);lEleL.push_back("PCB");
  lThickL.push_back(0.1*mm);lEleL.push_back("Si");
  lThickL.push_back(0.1*mm);lEleL.push_back("Si");
  lThickL.push_back(0.1*mm);lEleL.push_back("Si");

  std::vector<G4double> lThickR;
  std::vector<std::string> lEleR;
  // Rear module
  lThickR.push_back(0.6*mm*rScale);lEleR.push_back("WCu");
  lThickR.push_back(6*mm*rScale);lEleR.push_back("Cu");
  lThickR.push_back(0.6*mm*rScale);lEleR.push_back("WCu");
  lThickR.push_back(0.1*mm);lEleR.push_back("Si");
  lThickR.push_back(0.1*mm);lEleR.push_back("Si");
  lThickR.push_back(0.1*mm);lEleR.push_back("Si");
  lThickR.push_back(pcbThick);lEleR.push_back("PCB");
  lThickR.push_back(airThick);lEleR.push_back("Air");

  // Forward EE
  for(unsigned i=0; i<nFEModsNet; i++) {
    m_caloStruct.push_back( SamplingSection(lThickL,lEleL) );
    m_caloStruct.push_back( SamplingSection(lThickR,lEleR) );
  }

  // Middle EE
  unsigned nMEModsNom = 5;
  unsigned nMEModsSub = 0;
  unsigned nMEModsNet = nMEModsNom - nMEModsSub;
  rScale = ((double) nMEModsNom)/((double) nMEModsNet);
  lThickL[0] = 0.5*mm*rScale;
  lThickL[2] = 2.8*mm*rScale;//new
  lThickL[4] = 0.5*mm*rScale;
  lThickR[0] = 1.2*mm*rScale;//new
  lThickR[1] = 6.0*mm*rScale;
  lThickR[2] = 1.2*mm*rScale;//new
  for(unsigned i=0; i<nMEModsNet; i++) {
    m_caloStruct.push_back( SamplingSection(lThickL,lEleL) );
    m_caloStruct.push_back( SamplingSection(lThickR,lEleR) );
  }

  // Back EE
  unsigned nBEModsNom = 4;
  unsigned nBEModsSub = 1;
  unsigned nBEModsNet = nBEModsNom - nBEModsSub;
  rScale = ((double) nBEModsNom)/((double) nBEModsNet);
  lThickL[0] = 0.5*mm*rScale;
  lThickL[2] = 4.2*mm*rScale;//new
  lThickL[4] = 0.5*mm*rScale;
  lThickR[0] = 2.2*mm*rScale;//new
  lThickR[1] = 6.0*mm*rScale;
  lThickR[2] = 2.2*mm*rScale;//new
  for(unsigned i=0; i<nBEModsNet; i++) {
    m_caloStruct.push_back( SamplingSection(lThickL,lEleL) );
    if (i==(nBEModsNet-1)) {lThickR.push_back(0.5*mm);lEleR.push_back("Cu");}
    m_caloStruct.push_back( SamplingSection(lThickR,lEleR) );
  }  
	
	if(version_==v_HGCAL_v6){
	  //add HCAL
	  buildHGCALFHE(6);
	  buildHGCALBHE(6);
	}

	break;
      }

    case v_HGCALEE_v5: case v_HGCALEE_v5_gap4: case v_HGCAL_v5: case v_HGCAL_v5_gap4:
      {
	G4cout << "[DetectorConstruction] starting v_HCALEE_v5"<< G4endl;
	G4double airThick = 2*mm;
	if(version_==v_HGCALEE_v5_gap4 || version_==v_HGCAL_v5_gap4) airThick = 4*mm;

	G4double pcbThick = 1.2*mm;

	//first and last layers
	std::vector<G4double> lThick1;
	std::vector<std::string> lEle1;
	lThick1.push_back(180.*mm);lEle1.push_back("NeutMod");
	lThick1.push_back(0*mm);lEle1.push_back("Cu");
	lThick1.push_back(0*mm);lEle1.push_back("W");
	lThick1.push_back(0.5*mm);lEle1.push_back("Cu");
	lThick1.push_back(airThick);lEle1.push_back("Air");
	lThick1.push_back(pcbThick);lEle1.push_back("PCB");
	lThick1.push_back(0.1*mm);lEle1.push_back("Si");
	lThick1.push_back(0.1*mm);lEle1.push_back("Si");
	lThick1.push_back(0.1*mm);lEle1.push_back("Si");
	m_caloStruct.push_back( SamplingSection(lThick1,lEle1) );

	std::vector<G4double> lThickL;
	std::vector<std::string> lEleL;
	std::vector<G4double> lThickR;
	std::vector<std::string> lEleR;
	for (unsigned i=0;i<14;i++){
	  if (dropLayer_[2*i]) {
	    lThickR.clear();lEleR.clear();
	    lThickR.push_back(3*mm);lEleR.push_back("Cu");
	    lThickR.push_back(1*mm);lEleR.push_back("Pb");
	    //reset to 0.5 Cu and no lead for following layers
	    if (i>0) {
	      lThickR[0] = 0.5;
	      lThickR[1] = 0;
	    }
	    lThickR.push_back(absThickW_[i]);lEleR.push_back("W");
	    lThickR.push_back(0.5*mm);lEleR.push_back("Cu");
	    lThickR.push_back(airThick);lEleR.push_back("Air");
	    lThickR.push_back(absThickPb_[i]);lEleR.push_back("Pb");
	    lThickR.push_back(3*mm);lEleR.push_back("Cu");
	    lThickR.push_back(0.1*mm);lEleR.push_back("Si");
	    lThickR.push_back(0.1*mm);lEleR.push_back("Si");
	    lThickR.push_back(0.1*mm);lEleR.push_back("Si");
	    lThickR.push_back(pcbThick);lEleR.push_back("PCB");
	    lThickR.push_back(airThick);lEleR.push_back("Air");
	    m_caloStruct.push_back( SamplingSection(lThickR,lEleR) );
	  }
	  else {
	    lThickL.push_back(3*mm);lEleL.push_back("Cu");
	    lThickL.push_back(1*mm);lEleL.push_back("Pb");
	    //reset to 0.5 Cu and no lead for following layers
	    if (i>0) {
	      if (lThickL.size()==2){
		lThickL[0] = 0.5;
		lThickL[1] = 0;
	      }
	      else {
		lThickL[3] = 0.5;
		lThickL[4] = 0;
	      }
	    }
	    lThickL.push_back(absThickW_[i]);lEleL.push_back("W");
	    lThickL.push_back(0.5*mm);lEleL.push_back("Cu");
	    lThickL.push_back(airThick);lEleL.push_back("Air");
	    lThickL.push_back(pcbThick);lEleL.push_back("PCB");
	    lThickL.push_back(0.1*mm);lEleL.push_back("Si");
	    lThickL.push_back(0.1*mm);lEleL.push_back("Si");
	    lThickL.push_back(0.1*mm);lEleL.push_back("Si");
	    m_caloStruct.push_back( SamplingSection(lThickL,lEleL) );
	    lThickL.clear();lEleL.clear();
	    if (dropLayer_[2*i+1]) {
	      lThickL.clear();lEleL.clear();
	      lThickL.push_back(3*mm);lEleL.push_back("Cu");
	      lThickL.push_back(absThickPb_[i]);lEleL.push_back("Pb");
	      lThickL.push_back(airThick);lEleL.push_back("Air");
	    }
	    else {
	      lThickR.clear();lEleR.clear();
	      lThickR.push_back(3*mm);lEleR.push_back("Cu");
	      lThickR.push_back(absThickPb_[i]);lEleR.push_back("Pb");
	      lThickR.push_back(3*mm);lEleR.push_back("Cu");
	      lThickR.push_back(0.1*mm);lEleR.push_back("Si");
	      lThickR.push_back(0.1*mm);lEleR.push_back("Si");
	      lThickR.push_back(0.1*mm);lEleR.push_back("Si");
	      lThickR.push_back(pcbThick);lEleR.push_back("PCB");
	      lThickR.push_back(airThick);lEleR.push_back("Air");
	      m_caloStruct.push_back( SamplingSection(lThickR,lEleR) );
	    }
	  }
	}

	/*	//lThickR.push_back(0.5*mm);lEleR.push_back("Cu");
	for(unsigned i=0; i<4; i++) {
	  m_caloStruct.push_back( SamplingSection(lThickL,lEleL) );
	  m_caloStruct.push_back( SamplingSection(lThickR,lEleR) );
	}
	lThickL[2] = 2.8*mm;
	lThickR[1] = 2.1*mm;
	for(unsigned i=0; i<5; i++) {
	  m_caloStruct.push_back( SamplingSection(lThickL,lEleL) );
	  m_caloStruct.push_back( SamplingSection(lThickR,lEleR) );
	}
	lThickL[2] = 4.2*mm;
	lThickR[1] = 4.4*mm;
	for(unsigned i=0; i<4; i++) {
	  m_caloStruct.push_back( SamplingSection(lThickL,lEleL) );
	  m_caloStruct.push_back( SamplingSection(lThickR,lEleR) );
	}
	*/

	//last layer: add Cu+W+Cu...
	lThick1[0] = 0.*mm;
	lThick1[1] = 0.5*mm;
	lThick1[2] = absThickW_[14];
	//add last structure layers
	lThick1.push_back(3*mm);lEle1.push_back("Cu");
	lThick1.push_back(1*mm);lEle1.push_back("Pb");
	m_caloStruct.push_back( SamplingSection(lThick1,lEle1) );

	if(version_==v_HGCAL_v5 || version_==v_HGCAL_v5_gap4){
	  //add HCAL
	  buildHGCALFHE(5);
	  buildHGCALBHE(5);
	}

	break;
      }


    case v_HGCALEE_Si80: case v_HGCALEE_Si120: case v_HGCALEE_Si200: case v_HGCALEE_Si500: case v_HGCALEE_gap1: case  v_HGCALEE_CALICE: case v_HGCALEE_inverted: case v_HGCALEE_concept: case v_HGCALEE_W: case v_HGCALEE_gap4: case v_HGCALEE_prePCB: case v_HGCAL:
      {
	float siWidth(0.300), gap(2),pad(2);
	float pb1(1.60), pb2(3.30), pb3(5.60), cu(3.0);
	unsigned n1(10),n2(10),n3(10);

	if(version_==v_HGCALEE_Si80)     {siWidth=0.08;                               }
	if(version_==v_HGCALEE_Si120)    {siWidth=0.120;                              }
	if(version_==v_HGCALEE_Si500)    {siWidth=0.500;                              }
	if(version_==v_HGCALEE_gap1)     {gap=1;                                      }
	if(version_==v_HGCALEE_gap4)     {gap=4;                                      }
	if(version_==v_HGCALEE_CALICE)   {pb1=1.07;                                   }
	if(version_==v_HGCALEE_inverted) {pb2=1.63; pb1=3.32;                         }
	if(version_==v_HGCALEE_W)        {pb1=0.70;  pb2=2.10; pb3=3.50;}
	if(version_==v_HGCALEE_prePCB)   {addPrePCB_=true;}
	G4cout << "[DetectorConstruction] starting v_HGCAL* with Si width=" << siWidth << " and gap " << gap << G4endl;
	
	//add 1um silicon to track incoming particles
	//if(version_==v_HGCALEE_SiDummy) m_caloStruct.push_back( SamplingSection(0,0,0.001*mm,0,0) );
	std::vector<G4double> lThick;
	std::vector<std::string> lEle;
	lThick.push_back(pb1);lEle.push_back(version_==v_HGCALEE_W?"W":"Pb");
	lThick.push_back(cu); lEle.push_back("Cu");
	//add PCB to shield from delta-rays?
	if (addPrePCB_) {
	  lThick.push_back(pad);lEle.push_back("PCB");
	}
	lThick.push_back(siWidth/3.);lEle.push_back("Si");
	lThick.push_back(siWidth/3.);lEle.push_back("Si");
	lThick.push_back(siWidth/3.);lEle.push_back("Si");
	lThick.push_back(pad);lEle.push_back("PCB");
	lThick.push_back(gap);lEle.push_back("Air");

	if(version_==v_HGCALEE_concept || version_==v_HGCAL) {
	  lThick[0] = 0;
	  //lThick[1] = 0;
	  m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	}

	lThick[0] = pb1;
	lThick[1] = cu;
	for(unsigned i=0; i<n1; i++)         m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	lThick[0] = pb2;
	for(unsigned i=0; i<n2; i++)         m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	lThick[0] = pb3;
	for(unsigned i=0; i<n3; i++)         m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	
	if (version_==v_HGCAL){
	  buildHGCALFHE(4);
	  buildHGCALBHE(4);
	}
	break;
      }
    case v_HGCALHE: case v_HGCALHE_CMSSWv4: case v_HGCALHE_v6:
      {
	//add HCAL
	if (version_== v_HGCALHE_CMSSWv4) {
	  buildHGCALFHE(41);
	  buildHGCALBHE(4);
	}
	else if (version_== v_HGCALHE_v6) {
	  buildHGCALFHE(6);
	  buildHGCALBHE(6);
	}
	else {
	  buildHGCALFHE(4);
	  buildHGCALBHE(4);
	}
	break;
      }
    case v_HGCALHEScint:
      {
	buildHGCALBHE(4);
	break;
      }
    case v_HGCALBE_v6:
      {
	buildHGCALBHE(6);
	break;
      }
    case v_HGCALHE_CALICE:
      {
	std::vector<G4double> lThick;
	std::vector<std::string> lEle;
	lThick.push_back(16.7*mm);lEle.push_back("Steel");
	lThick.push_back(1.25*mm);lEle.push_back("Air");
	lThick.push_back(2*mm);lEle.push_back("Steel");
	lThick.push_back(1.5*mm);lEle.push_back("CFMix");
	lThick.push_back(1*mm);lEle.push_back("PCB");
	lThick.push_back(0.115*mm);lEle.push_back("Polystyrole");
	lThick.push_back(5*mm);lEle.push_back("Scintillator");
	lThick.push_back(0.115*mm);lEle.push_back("Polystyrole");
	lThick.push_back(2*mm);lEle.push_back("Steel");
	lThick.push_back(1.25*mm);lEle.push_back("Air");
	//for(unsigned i=0; i<12; i++) m_caloStruct.push_back( SamplingSection(52,3.0,0.3,2.0,2.0) );
	//total 5.3 lambda = 22mm brass * 38 layers
	for(unsigned i=0; i<3; i++) m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	lThick[0] = 17.4*mm;
	for(unsigned i=0; i<21; i++) m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	lThick[0] = 17.6*mm;
	for(unsigned i=0; i<14; i++) m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	//last absorber chunk ??
	//m_caloStruct.push_back( SamplingSection(20.5,0.,0,0,0) );

	//for(unsigned i=0; i<38; i++) m_caloStruct.push_back( SamplingSection(21,0.,5,0.,0.) );
	lThick[0] = 21*mm;
	for(unsigned i=0; i<9; i++) m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	lThick[0] = 104*mm;
	for(unsigned i=0; i<7; i++) m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	break;
      }

    }

  DefineMaterials();
  SetMagField(0);
  m_detectorMessenger = new DetectorMessenger(this);
  UpdateCalorSize();
}

void DetectorConstruction::buildHGCALFHE(const unsigned aVersion){
  G4double airThick = 2*mm;
  if(version_==v_HGCAL_v5_gap4) airThick = 4*mm;
  std::vector<G4double> lThick;
  std::vector<std::string> lEle;
  if(aVersion==6) {
    airThick = 3*mm;
    G4double pcbthick = 2*mm;
    //putting all absorber in front of each Si layer to have correct reweighting 
    lThick.push_back(15.*mm);lEle.push_back("SSteel");
    lThick.push_back(45.*mm);lEle.push_back("Brass");
    lThick.push_back(0.5*mm); lEle.push_back("Cu");
    lThick.push_back(airThick);lEle.push_back("Air");
    lThick.push_back(pcbthick);lEle.push_back("PCB");
    lThick.push_back(0.1*mm);lEle.push_back("Si");
    lThick.push_back(0.1*mm);lEle.push_back("Si");
    lThick.push_back(0.1*mm);lEle.push_back("Si");
    m_caloStruct.push_back( SamplingSection(lThick,lEle) );

    lThick.clear();
    lEle.clear();
    lThick.push_back(1.);lEle.push_back("CFMix");
    lThick.push_back(6.); lEle.push_back("Cu");
    lThick.push_back(45.*mm);lEle.push_back("Brass");
    lThick.push_back(0.5*mm); lEle.push_back("Cu");
    lThick.push_back(airThick);lEle.push_back("Air");
    lThick.push_back(pcbthick);lEle.push_back("PCB");
    lThick.push_back(0.1*mm);lEle.push_back("Si");
    lThick.push_back(0.1*mm);lEle.push_back("Si");
    lThick.push_back(0.1*mm);lEle.push_back("Si");
    for(unsigned i=0; i<11; i++) {
      m_caloStruct.push_back( SamplingSection(lThick,lEle) );
    }
  }
  else {
    lThick.push_back(50.*mm);lEle.push_back("SSteel");
    if (aVersion==41) {lThick.push_back(52.*mm);lEle.push_back("Pb");}
    else {lThick.push_back(26.*mm);lEle.push_back("Brass");}
    lThick.push_back(3*mm); lEle.push_back("Cu");
    lThick.push_back(0.1*mm);lEle.push_back("Si");
    lThick.push_back(0.1*mm);lEle.push_back("Si");
    lThick.push_back(0.1*mm);lEle.push_back("Si");
    G4double pcbthick = (aVersion==4)? 2*mm : 1.2*mm;
    lThick.push_back(pcbthick);lEle.push_back("PCB");
    lThick.push_back(airThick);lEle.push_back("Air");
    
    for(unsigned i=0; i<12; i++) {
      //add an intermediate Si layer to study resolution improvement
      if (i>0) lThick[0] = 0;
      if (aVersion!=41){
	lThick[2] = 0;
	lThick[6] = 0;
	lThick[7] = 0;
	m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	lThick[0] = 0;
	lThick[2] = 3*mm;
	lThick[6] = pcbthick;
	lThick[7] = airThick;
      }
      if (i==11) lThick[7] = 103.6*mm;
      m_caloStruct.push_back( SamplingSection(lThick,lEle) );
    }
  }
}
//
void DetectorConstruction::buildHGCALBHE(const unsigned aVersion){
  std::vector<G4double> lThick;
  std::vector<std::string> lEle;
  //first layer
  lThick.push_back(1.*mm);lEle.push_back("CFMix");
  lThick.push_back(6.*mm); lEle.push_back("Cu");
  lThick.push_back(65.*mm);lEle.push_back("Air");
  lThick.push_back(2.*mm);lEle.push_back("Al");
  lThick.push_back(26.*mm);lEle.push_back("Foam");
  lThick.push_back(2.*mm);lEle.push_back("Al");
  lThick.push_back(78.*mm);lEle.push_back("Brass");
  if (aVersion==6) {
    lThick.push_back(2.6*mm);lEle.push_back("Air");
    lThick.push_back(3.8*mm);lEle.push_back("Scintillator");
    lThick.push_back(2.6*mm);lEle.push_back("Air");
  }
  else {
    lThick.push_back(9.*mm);lEle.push_back("Scintillator");
  }
  m_caloStruct.push_back( SamplingSection(lThick,lEle) );

  //other layers
  lThick.clear();
  lEle.clear();
  lThick.push_back(78.*mm);lEle.push_back("Brass");
  if (aVersion==6) {
    lThick.push_back(2.6*mm);lEle.push_back("Air");
    lThick.push_back(3.8*mm);lEle.push_back("Scintillator");
    lThick.push_back(2.6*mm);lEle.push_back("Air");
  }
  else {
    lThick.push_back(9.*mm);lEle.push_back("Scintillator");
  }

  unsigned maxi = (aVersion==4)?9:11;
  for(unsigned i=0; i<maxi; i++) {
    m_caloStruct.push_back( SamplingSection(lThick,lEle) );
  }
}
//
DetectorConstruction::~DetectorConstruction() { delete m_detectorMessenger;}

//
void DetectorConstruction::DefineMaterials()
{ 
  G4NistManager* nistManager = G4NistManager::Instance();
  m_materials["Abs"] = (version_== v_CALICE || version_==v_HGCALEE_W) ? 
    nistManager->FindOrBuildMaterial("G4_W",false) :
    nistManager->FindOrBuildMaterial("G4_Pb",false);
  m_materials["Al"] = nistManager->FindOrBuildMaterial("G4_Al",false); 
  m_materials["W"] = nistManager->FindOrBuildMaterial("G4_W",false); 
  m_materials["Pb"] = nistManager->FindOrBuildMaterial("G4_Pb",false); 
  m_materials["Cu"] = nistManager->FindOrBuildMaterial("G4_Cu",false); 
  m_materials["Si"] = nistManager->FindOrBuildMaterial("G4_Si",false);
  m_materials["Zn"] = nistManager->FindOrBuildMaterial("G4_Zn",false);
  m_materials["Air"]=nistManager->FindOrBuildMaterial("G4_AIR",false);
  m_materials["Fe"] = nistManager->FindOrBuildMaterial("G4_Fe",false);
  m_materials["Mn"] = nistManager->FindOrBuildMaterial("G4_Mn",false);
  m_materials["C"] = nistManager->FindOrBuildMaterial("G4_C",false); 
  m_materials["H"] = nistManager->FindOrBuildMaterial("G4_H",false); 
  m_materials["Cl"] = nistManager->FindOrBuildMaterial("G4_Cl",false); 
  m_materials["Cr"] = nistManager->FindOrBuildMaterial("G4_Cr",false); 
  m_materials["Ni"] = nistManager->FindOrBuildMaterial("G4_Ni",false); 

  m_materials["PCB"] = new G4Material("G10",1.700*g/cm3,4);
  m_materials["PCB"]->AddElement(nistManager->FindOrBuildElement(14), 1);
  m_materials["PCB"]->AddElement(nistManager->FindOrBuildElement(8) , 2);
  m_materials["PCB"]->AddElement(nistManager->FindOrBuildElement(6) , 3);
  m_materials["PCB"]->AddElement(nistManager->FindOrBuildElement(1) , 3);
  m_materials["Brass"]= new G4Material("Brass",8.5*g/cm3,2);
  m_materials["Brass"]->AddMaterial(m_materials["Cu"]  , 70*perCent);
  m_materials["Brass"]->AddMaterial(m_materials["Zn"]  , 30*perCent);
  m_materials["Steel"]= new G4Material("Steel",7.87*g/cm3,3);
  m_materials["Steel"]->AddMaterial(m_materials["Fe"]  , 0.9843);
  m_materials["Steel"]->AddMaterial(m_materials["Mn"], 0.014);
  m_materials["Steel"]->AddMaterial(m_materials["C"], 0.0017);
  m_materials["SSteel"]= new G4Material("SSteel",8.02*g/cm3,4);
  m_materials["SSteel"]->AddMaterial(m_materials["Fe"]  , 0.70);
  m_materials["SSteel"]->AddMaterial(m_materials["Mn"], 0.01);
  m_materials["SSteel"]->AddMaterial(m_materials["Cr"], 0.19);
  m_materials["SSteel"]->AddMaterial(m_materials["Ni"], 0.10);
  m_materials["AbsHCAL"] = (version_== v_HGCALHE_CALICE) ?
    m_materials["Steel"]:
    m_materials["Brass"];
  m_materials["Scintillator"]= nistManager->FindOrBuildMaterial("G4_POLYSTYRENE",false); 
  //m_materials["Scintillator"]= new G4Material("Scintillator",1.032*g/cm3,2);
  //m_materials["Scintillator"]->AddMaterial(m_materials["C"]  , 91.512109*perCent);
  //m_materials["Scintillator"]->AddMaterial(m_materials["H"]  , 8.4878906*perCent);
  G4cout << m_materials["Scintillator"] << G4endl;
  m_materials["Polystyrole"]= new G4Material("Polystyrole",1.065*g/cm3,2);
  m_materials["Polystyrole"]->AddMaterial(m_materials["H"]  , 50*perCent);
  m_materials["Polystyrole"]->AddMaterial(m_materials["C"]  , 50*perCent);

  m_materials["PVC"]= new G4Material("PVC",1.350*g/cm3,3);
  m_materials["PVC"]->AddMaterial(m_materials["H"]  , 50*perCent);
  m_materials["PVC"]->AddMaterial(m_materials["C"]  , 33.33*perCent);
  m_materials["PVC"]->AddMaterial(m_materials["Cl"]  , 16.67*perCent);

  m_materials["CFMix"]= new G4Material("CFMix",0.120*g/cm3,3);
  m_materials["CFMix"]->AddMaterial(m_materials["Air"]  , 0.009);
  m_materials["CFMix"]->AddMaterial(m_materials["PVC"]  , 0.872);
  m_materials["CFMix"]->AddMaterial(m_materials["Polystyrole"]  , 0.119);

  m_materials["Foam"]= new G4Material("Foam",0.0999*g/cm3,2);
  m_materials["Foam"]->AddMaterial(m_materials["C"]  , 0.856);
  m_materials["Foam"]->AddMaterial(m_materials["H"]  , 0.144);

  m_materials["WCu"]= new G4Material("WCu",14.979*g/cm3,2);
  m_materials["WCu"]->AddMaterial(m_materials["W"]  , 75*perCent);
  m_materials["WCu"]->AddMaterial(m_materials["Cu"]  , 25*perCent);

  m_materials["NeutMod"]= new G4Material("NeutMod",0.950*g/cm3,2);
  m_materials["NeutMod"]->AddMaterial(m_materials["C"]  , 0.85628);
  m_materials["NeutMod"]->AddMaterial(m_materials["H"]  , 0.14372);

}

//
void DetectorConstruction::UpdateCalorSize(){  

  m_CalorSizeZ=0;
  for(size_t i=0; i<m_caloStruct.size(); i++)
    m_CalorSizeZ=m_CalorSizeZ+m_caloStruct[i].Total_thick;

  if (model_ == DetectorConstruction::m_SIMPLE_20)
    m_CalorSizeXY=200;
  else if (model_ == DetectorConstruction::m_SIMPLE_50)
    m_CalorSizeXY=500;
  else if (model_ == DetectorConstruction::m_SIMPLE_100)
    m_CalorSizeXY=1000;
  else if (model_ == DetectorConstruction::m_FULLSECTION){
    m_CalorSizeXY=2800;//1700;
    m_minRadius = 150;
    m_maxRadius = m_CalorSizeXY;
    m_minEta = 1.4;
    m_maxEta = 3.7;
    m_z0pos = 2990;//3170;
  }
  else m_CalorSizeXY=200;

  m_WorldSizeZ=m_CalorSizeZ*1.1;  
  m_WorldSizeXY=m_CalorSizeXY*1.1;

  if (model_ == DetectorConstruction::m_FULLSECTION) 
    G4cout << "[DetectorConstruction][UpdateCalorSize] Z x minR * maxR = " 
	   << m_CalorSizeZ << " x " 
	   << m_minRadius << " x " 
	   << m_maxRadius 
	   << " mm, eta range "
	   << m_minEta << " - " 
	   << m_maxEta
	   << G4endl;
  else G4cout << "[DetectorConstruction][UpdateCalorSize] Z x XY = " 
	      << m_CalorSizeZ << " x " 
	      << m_CalorSizeXY << " mm " <<  G4endl;

}

//
G4VPhysicalVolume* DetectorConstruction::Construct()
{

  //clean old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  //world
  G4double expHall_z = 6*m;
  G4double expHall_x = 3*m;
  G4double expHall_y = 3*m;

  G4Box* experimentalHall_box = new G4Box("expHall_box",expHall_x,expHall_y,expHall_z);

  G4LogicalVolume* experimentalHall_log = new G4LogicalVolume(experimentalHall_box, m_materials["Air"],"expHall_log");
  G4VPhysicalVolume* experimentalHall_phys
    = new G4PVPlacement(0,                       // no rotation
			G4ThreeVector(0.,0.,0.), // translation position
			experimentalHall_log,    // its logical volume
			"expHall",               // its name
			0,                       // its mother volume
			false,                   // no boolean operations
			0);                      // its copy number

  //detector's World
  G4double pos_x = 0.;
  G4double pos_y = 0.;
  G4double pos_z = 0.;
  if (model_ == DetectorConstruction::m_FULLSECTION){
    pos_x = 0.;
    pos_y = 0.;
    pos_z = m_z0pos+m_CalorSizeZ/2;
  }

  if (model_ == DetectorConstruction::m_FULLSECTION){
    m_solidWorld = new G4Tubs("Wbox",m_minRadius*1.1,m_maxRadius*1.1,m_WorldSizeZ/2,0,2*pi);
  }
  else {
    m_solidWorld = new G4Box("Wbox",m_WorldSizeXY/2,m_WorldSizeXY/2,m_WorldSizeZ/2);
  }
  m_logicWorld = new G4LogicalVolume(m_solidWorld, m_materials["Air"], "Wlog");
  m_physWorld = new G4PVPlacement(0, G4ThreeVector(pos_x,pos_y,pos_z), m_logicWorld, "Wphys", experimentalHall_log, false, 0);



  //build the stack
  G4double zOffset(-m_CalorSizeZ/2), zOverburden(0.);
  char nameBuf[10];
  G4CSGSolid *solid;

  G4double totalLengthX0 = 0;
  G4double totalLengthL0 = 0;

  for(size_t i=0; i<m_caloStruct.size(); i++)
    {

      const unsigned nEle = m_caloStruct[i].n_elements;
      //index for counting Si sensitive layers
      unsigned idx = 0;

      for (unsigned ie(0); ie<nEle;++ie){
	std::string eleName = m_caloStruct[i].ele_name[ie];
	sprintf(nameBuf,"%s%d",eleName.c_str(),int(i+1));
	if (eleName=="Si") {
	  sprintf(nameBuf,"Si%d_%d",int(i+1),idx); 
	  idx++;
	}
	std::string baseName(nameBuf);
	G4double thick = m_caloStruct[i].ele_thick[ie];
	if(thick>0){
	  solid = constructSolid(baseName,thick,zOffset+zOverburden);
	  G4LogicalVolume *logi = new G4LogicalVolume(solid, m_materials[eleName], baseName+"log");
	  m_caloStruct[i].ele_X0[ie] = m_materials[eleName]->GetRadlen();
	  m_caloStruct[i].ele_L0[ie] = m_materials[eleName]->GetNuclearInterLength();
	  G4cout << "************ " << eleName << " layer " << i << " X0=" << m_caloStruct[i].ele_X0[ie] << " L0=" << m_caloStruct[i].ele_L0[ie] << " w=" << m_caloStruct[i].ele_thick[ie] << "mm, d=" << m_materials[eleName]->GetDensity();
	  //G4cout << G4endl;
	  //G4cout << *(m_materials[eleName]->GetMaterialTable()) << G4endl;

	  totalLengthX0 += m_caloStruct[i].ele_thick[ie]/m_caloStruct[i].ele_X0[ie];
	  G4cout << " TotX0=" << totalLengthX0;// << G4endl;
	  totalLengthL0 += m_caloStruct[i].ele_thick[ie]/m_caloStruct[i].ele_L0[ie];
	  G4cout << " TotLambda=" << totalLengthL0 << G4endl;

	  if (m_caloStruct[i].isSensitiveElement(ie)) m_logicSi.push_back(logi);
	  
	  m_caloStruct[i].ele_vol[ie]= new G4PVPlacement(0, G4ThreeVector(0.,0.,zOffset+zOverburden+thick/2), logi, baseName+"phys", m_logicWorld, false, 0);


	  G4VisAttributes *simpleBoxVisAtt= new G4VisAttributes(m_caloStruct[i].g4Colour(ie));
	  simpleBoxVisAtt->SetVisibility(true);
	  logi->SetVisAttributes(simpleBoxVisAtt);
	  zOverburden = zOverburden + thick;
	  //for sensitive volumes
	  //add region to be able to set specific cuts for it
	  //just for Si
	  if (eleName=="Si"){
	    unsigned nlogicsi = m_logicSi.size();
	    G4Region* aRegion = new G4Region(baseName+"Reg");
	    m_logicSi[nlogicsi-1]->SetRegion(aRegion);
	    aRegion->AddRootLogicalVolume(m_logicSi[nlogicsi-1]);
	  }
	}

      }//loop on elements
    }//loop on layers
      //                                        
  // Visualization attributes
  //
  m_logicWorld->SetVisAttributes (G4VisAttributes::Invisible);

  //return m_physWorld;
  return experimentalHall_phys;
}

//
void DetectorConstruction::SetMagField(G4double fieldValue)
{

  if(fieldValue<=0) return; 

  //apply a global uniform magnetic field along Z axis
  G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  if(m_magField) delete m_magField;                //delete the existing magn field
  m_magField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));
  fieldMgr->SetDetectorField(m_magField);
  fieldMgr->CreateChordFinder(m_magField);
  fieldMgr->SetDetectorField(m_magField);  
}

void DetectorConstruction::SetDetModel(G4int model)
{
  if (model <= 0) return;
  std::cout << " -- Setting detector model to " << model << std::endl;
  model_ = model;
}

void DetectorConstruction::SetWThick(std::string thick)
{
  if (thick.size() <= 0) return;
  std::cout << " -- Setting W thick to " << thick << std::endl;
  std::vector<std::string> vec;
  boost::split(vec, thick, boost::is_any_of(","));
  absThickW_.resize(vec.size(),0);
  for (unsigned iE(0); iE<vec.size(); ++iE){//loop on elements
    std::istringstream(vec[iE])>>absThickW_[iE];
    std::cout << absThickW_[iE] << " ";
  }
  std::cout << std::endl;
}

void DetectorConstruction::SetPbThick(std::string thick)
{
  if (thick.size() <= 0) return;
  std::cout << " -- Setting Pb thick to " << thick << std::endl;
  std::vector<std::string> vec;
  boost::split(vec, thick, boost::is_any_of(","));
  absThickPb_.resize(vec.size(),0);
  for (unsigned iE(0); iE<vec.size(); ++iE){//loop on elements
    std::istringstream(vec[iE])>>absThickPb_[iE];
    std::cout << absThickPb_[iE] << " ";
  }
  std::cout << std::endl;
}

void DetectorConstruction::SetDropLayers(std::string layers)
{
  dropLayer_.resize(28,false);
  if (layers.size() <= 0) return;
  std::cout << " -- Dropping layers " << layers << std::endl;
  std::vector<std::string> vec;
  boost::split(vec, layers, boost::is_any_of(","));
  for (unsigned iE(0); iE<vec.size(); ++iE){//loop on elements
    unsigned layerId = 0;
    std::istringstream(vec[iE])>>layerId;
    if (layerId>0 && layerId<29) dropLayer_[layerId-1] = true;
    else std::cout << " -- invalid layer to drop, ignoring..." << std::endl;
  }
  for (unsigned iE(0); iE<dropLayer_.size(); ++iE){//loop on elements
    std::cout << dropLayer_[iE] << " ";
  }
  std::cout << std::endl;
}

G4CSGSolid *DetectorConstruction::constructSolid (std::string baseName, G4double thick, G4double zpos){
  
  G4CSGSolid *solid;
  if (model_ == DetectorConstruction::m_FULLSECTION){
    double minR = tan(2*atan(exp(-m_maxEta)))*(zpos+m_z0pos+m_CalorSizeZ/2);
    double maxR = tan(2*atan(exp(-m_minEta)))*(zpos+m_z0pos+m_CalorSizeZ/2);
    //std::cout << " zpos = " << zpos+m_z0pos+m_CalorSizeZ/2 << " radius range " << minR << " " << maxR << std::endl;
    solid = new G4Tubs(baseName+"box",minR,maxR,thick/2,0,2*pi);
  }
  else{
    solid = new G4Box(baseName+"box", m_CalorSizeXY/2, m_CalorSizeXY/2, thick/2 );
  }
  return solid;
}
