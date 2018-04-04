#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "SamplingSection.hh"

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include <map>
#include <string>

#include "TH1F.h"

class G4VSolid;
class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class DetectorMessenger;
class G4Colour;

/**
   @class DetectorConstruction
   @short builds a simple detector
 */
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

  enum DetectorVersion { 
    v_CALICE=0,
    v_HGCALEE_Si80=1,
    v_HGCALEE_Si120=2,
    v_HGCALEE_Si200=3,
    v_HGCALEE_Si500=4,
    v_HGCALEE_gap1=5,
    v_HGCALEE_CALICE=6,
    v_HGCALEE_inverted=7,
    v_HGCALEE_concept=8,
    v_HGCALEE_W=9,
    v_HGCALEE_gap4=10,
    v_HGCALEE_prePCB=11,
    v_HGCALEE_v5=12,
    v_HGCALEE_v5_gap4=13,
    v_HGCAL=20,
    v_HGCALHE=21,
    v_HGCALHEScint=22,
    v_HGCALHE_CALICE=23,
    v_HGCALHE_CMSSWv4=24,
    v_HGCAL_v5=25,
    v_HGCAL_v5_gap4=26,
    v_HGCALHE_v5=27,
    v_HGCALBE_v5=28,
    v_HGCALEE_v6=30,
    v_HGCALHE_v6=31,
    v_HGCALBE_v6=32,
    v_HGCAL_v6=33,
    v_HGCALEE_v624=34,
    v_HGCALEE_v618=35,
    v_HGCAL_v624=36,
    v_HGCAL_v618=37,
    v_HGCALHE_v624=38,
    v_HGCALHE_v618=39,
    v_HGCALEE_v7=40,
    v_HGCALHE_v7=41,
    v_HGCALBE_v7=42,
    v_HGCAL_v7=43,
    v_HGCALHF=50,
    v_HGCAL_v7_HF=51,
    v_HGCALCu_HF=52,
    v_HGCALEE_v8=60,
    v_HGCALHE_v8=61,
    v_HGCALBE_v8=62,
    v_HGCAL_v8     = 63,  //TDR
    v_HGCALSci_v8  = 630, //TDR EE+FH Si air
    v_HGCAL_v9     = 64,  //24 uniform
    v_HGCALSci_v9  = 640, //24 Si air
    v_HGCAL_v10    = 65,  //9+9
    v_HGCALSci_v10 = 650, //9+9 Si air
    v_HGCAL_v11    = 66,  //18 uniform
    v_HGCALSci_v11 = 660, //18 Si air
    v_HGCALEE_TB=100,
    v_HGCALEE_TB_gap0=101,
    v_HGCALEE_TB_allW=102,
    v_HGCALEE_TB_samedEdx=103,
    v_HGCAL_2016TB=110
  };

  enum DetectorModel {
    m_SIMPLE_20=0,
    m_SIMPLE_50=1,
    m_FULLSECTION=2,
    m_SIMPLE_100=3,
    m_BOXWITHCRACK_100=4,
    m_2016TB=5
  };

  /**
     @short CTOR
   */
  DetectorConstruction(G4int ver=DetectorConstruction::v_CALICE, 
		       G4int mod=DetectorConstruction::m_SIMPLE_20,
		       G4int shape=1,
		       std::string absThickW="1.75,1.75,1.75,1.75,1.75,2.8,2.8,2.8,2.8,2.8,4.2,4.2,4.2,4.2,4.2",
		       std::string absThickPb="1,1,1,1,1,2.1,2.1,2.1,2.1,2.1,4.4,4.4,4.4,4.4",
		       std::string dropLayer="");

  void buildHGCALFHE(const unsigned aVersion,bool doAir=false);
  void buildHGCALBHE(const unsigned aVersion);
  void buildHF();
  /**
     @short calorimeter structure (sampling sections)
   */
  std::vector<SamplingSection> m_caloStruct;
  std::vector<SamplingSection> *getStructure() { return &m_caloStruct; }

  int getModel() const { return model_; }
  int getVersion() const { return version_; }
  unsigned getShape() const { return shape_; }



  const std::vector<G4LogicalVolume*>  & getSiLogVol() {return m_logicSi; }
  const std::vector<G4LogicalVolume*>  & getAlLogVol() {return m_logicAl; }
  const std::vector<G4LogicalVolume*>  & getAbsLogVol() {return m_logicAbs; }


  /**
     @short define the calorimeter materials
   */
  void DefineMaterials(); 
  std::map<std::string, G4Material *> m_materials;
  std::map<std::string, G4double > m_dEdx;
  std::map<std::string, G4Colour > m_colours;

  /**
     @short set magnetic field
   */
  void SetMagField(G4double fieldValue);
  G4UniformMagField* m_magField;      //pointer to the magnetic field

  /**
     @short set detector model
   */

  void SetDetModel(G4int model);

  void SetWThick(std::string thick);
  void SetPbThick(std::string thick);
  void SetDropLayers(std::string layers);

  /**
     @short DTOR
   */
  ~DetectorConstruction();
  

  /**
     @short getters
   */
  G4double GetCalorSizeXY() { return m_CalorSizeXY; }
  G4double GetCalorSizeZ()  { return m_CalorSizeZ; }
  G4double GetWorldSizeXY() { return m_WorldSizeXY; }
  G4double GetWorldSizeZ()  { return m_WorldSizeZ; }

  G4double GetMinEta() { return m_minEta0; }
  G4double GetMaxEta() { return m_maxEta0; }

  TH1F *getSensitiveZHisto(){ return sensitiveZ_h_; }
  TH1F *getEtaBoundaryHisto(bool isMin) { return isMin ? etaBoundaryMin_h_ : etaBoundaryMax_h_; }

  /**
     @short build the detector
   */

  G4VPhysicalVolume* Construct();

private:

  //detector version
  int version_;
  //integer to define detector model
  int model_;
  //integer to define shape of simhits: 1=hexa, 2=diamonds, 3=triangles, 4=squares
  unsigned shape_;

  //add a pre PCB plate
  bool addPrePCB_;

  bool doHF_;
  unsigned firstHFlayer_;
  unsigned firstMixedlayer_;
  unsigned firstScintlayer_;
  unsigned nbhLayers_;

  std::vector<G4double> absThickW_;
  std::vector<G4double> absThickPb_;
  std::vector<G4bool> dropLayer_;

  double getEtaFromRZ(const double & r, const double & z);

  /**
     @short compute the calor dimensions
   */
  void UpdateCalorSize(); 

  /**
     @short build the calorimeter
   */
  G4VPhysicalVolume* ConstructCalorimeter();     

  void buildSectorStack(const unsigned sectorNum,
			const G4double & minL, 
			const G4double & width);

  void fillInterSectorSpace(const unsigned sectorNum,
			    const G4double & minL, 
			    const G4double & width);

  G4double getCrackOffset(size_t layer);
  G4double getAngOffset(size_t layer);

  G4VSolid *constructSolid (const unsigned layer, std::string baseName, G4double thick, G4double zpos,const G4double & minL, const G4double & width, const bool isHF=false);
  G4VSolid *constructSolid (std::string baseName, G4double thick, G4double zpos,const G4double & minL, const G4double & width, const double & etamin, const double & etamax);

  G4VSolid *constructSupportCone (std::string baseName, G4double thick, G4double zpos,const G4double & minL, const G4double & width, const double & etamin, const double & etamax);
  G4VSolid *constructSupportCone (std::string baseName, G4double thick, G4double zpos,const G4double & minL, const G4double & width);

  std::vector<G4Material* > m_SensitiveMaterial;
  
  G4double           m_CalorSizeXY, m_CalorSizeZ;
  G4double           m_minRadius,m_maxRadius;
  G4double           m_minRadiusHF,m_maxRadiusHF;
  //define per layer
  std::vector<G4double>           m_minEta,m_maxEta;
  G4double           m_minEta0,m_maxEta0;
  G4double           m_minEtaHF,m_maxEtaHF;
  G4double           m_z0pos, m_z0HF;
  G4double           m_WorldSizeXY, m_WorldSizeZ;
  G4double m_nSectors,m_sectorWidth,m_interSectorWidth;
            
  G4VSolid*          m_solidWorld;    //pointer to the solid World 
  G4LogicalVolume*   m_logicWorld;    //pointer to the logical World
  G4VPhysicalVolume* m_physWorld;     //pointer to the physical World  
  
  std::vector<G4LogicalVolume*>   m_logicSi;    //pointer to the logical Si volumes
  std::vector<G4LogicalVolume*>   m_logicAl;    //pointer to the logical Si volumes
  std::vector<G4LogicalVolume*>   m_logicAbs;    //pointer to the logical absorber volumes situated just before the si

  DetectorMessenger* m_detectorMessenger;  //pointer to the Messenger

  TH1F *sensitiveZ_h_,*etaBoundaryMin_h_,*etaBoundaryMax_h_;

};


#endif

