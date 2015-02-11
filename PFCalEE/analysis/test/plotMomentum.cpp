#include <iostream>
#include <sstream>
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"

#include "TH1F.h"
#include "TFile.h"

using namespace HepMC;

int main(int argc, char ** argv) {

  std::cout << " -- Starting program..." << std::endl;

  TFile *fout = TFile::Open("PLOTS/output_hepmc.root","RECREATE");
  fout->cd();
  TH1F *hProton = new TH1F("hProton",";p (GeV);particles",1000,0,100);
  TH1F *hNeutron = new TH1F("hNeutron",";p (GeV);particles",1000,0,100);
  TH1F *hPipm = new TH1F("hPipm",";p (GeV);particles",1000,0,100);
  TH1F *hProtonLog = new TH1F("hProtonLog",";log(p) (log(GeV));particles",100,-2,2);
  TH1F *hNeutronLog = new TH1F("hNeutronLog",";log(p) (log(GeV));particles",100,-2,2);
  TH1F *hPipmLog = new TH1F("hPipmLog",";log(p) (log(GeV));particles",100,-2,2);

  const unsigned nFiles = 6;

  for (unsigned i(0);i<nFiles;++i){

    // specify an input file
    std::ostringstream lname;
    lname << "/afs/cern.ch/work/p/pdauncey/public/Pythia140305_";
    if (i==0) lname << "000000";
    else if (i<10) lname << "00000" << i;
    else if (i<100) lname << "0000" << i;
    else if (i<1000) lname << "000" << i;
    else if (i<10000) lname << "00" << i;
    else if (i<100000) lname << "0" << i;
    else lname << i;
    lname << ".dat";
    HepMC::IO_GenEvent ascii_in(lname.str().c_str(),std::ios::in);
    
    // get the first event
    HepMC::GenEvent* evt = ascii_in.read_next_event();
    // loop until we run out of events
    while ( evt ) {
      unsigned ievt =  evt->event_number();
      if (ievt%100==0) std::cout << "Processing Event Number "
				 << ievt
				 << std::endl;
      // analyze the event
      HepMC::GenEvent::particle_const_iterator lPart = evt->particles_begin();
      //unsigned counter = 0;
      
      //std::cout << " -- Number of particles: " << evt->particles_size()
      //<< std::endl;
      
      for (; lPart!=evt->particles_end();++lPart){
	//std::cout << counter << " " 
	//<< (*lPart)->pdg_id() <<  " " 
	//<< (*lPart)->status() 
	//<< std::endl;
	//counter++;
	if ((*lPart)->status()!=1) continue;
	double p = sqrt(pow((*lPart)->momentum().px(),2)+pow((*lPart)->momentum().py(),2)+pow((*lPart)->momentum().pz(),2));
	//if (fabs((*lPart)->momentum().eta())>2.8 && fabs((*lPart)->momentum().eta())<3.0){
	if (fabs((*lPart)->momentum().eta())<2.5){
	  if (fabs((*lPart)->pdg_id())==2212) {hProton->Fill(p);hProtonLog->Fill(log10(p));}
	  else if (fabs((*lPart)->pdg_id())==2112) {hNeutron->Fill(p);hNeutronLog->Fill(log10(p));}
	  else if (fabs((*lPart)->pdg_id())==211) {hPipm->Fill(p);hPipmLog->Fill(log10(p));}
	}
      }
      
      // delete the created event from memory
      delete evt;
      // read the next event
      ascii_in >> evt;
      ievt++;
    }
  }

  fout->Write();
  
  return 0;

}//main




