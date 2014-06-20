#ifndef EventAction_h
#define EventAction_h 1

#include "SamplingSection.hh"

#include "G4ThreeVector.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"

#include "TFile.h"
//#include "TNtuple.h"
#include "TTree.h"
#include "SamplingSection.hh"
#include "G4SiHit.hh"
#include "HGCSSEvent.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSGenParticle.hh"

#include <vector>
#include <map>

class RunAction;
class EventActionMessenger;

class EventAction : public G4UserEventAction
{
public:
  EventAction();
  virtual ~EventAction();
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);

  void Detect(G4double edep, G4double stepl,G4double globalTime, G4int pdgId, 
	      G4VPhysicalVolume *volume, const G4ThreeVector & position,
	      G4int trackID, G4int parentID,
	      const HGCSSGenParticle & genPart);
  
  //void Detect(G4double edep, G4double stepl,G4double globalTime, G4int pdgId, G4VPhysicalVolume *volume,int iyiz);

  void SetPrintModulo(G4int    val)  {printModulo = val;};
  void Add( std::vector<SamplingSection> *newDetector ) { detector_=newDetector; }
  //Float_t GetCellSize() { return cellSize_; }

  std::string GetFirstVolumeName() const;

private:
  RunAction*  runAct;
  std::vector<SamplingSection> *detector_;
  G4int     evtNb_,printModulo;
  TFile *outF_;
  TTree *tree_;
  HGCSSEvent event_;
  HGCSSSamplingSectionVec ssvec_;
  HGCSSSimHitVec hitvec_;
  HGCSSGenParticleVec genvec_;
  EventActionMessenger*  eventMessenger;

};

#endif

    
