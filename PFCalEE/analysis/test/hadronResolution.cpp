#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "boost/function.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TChain.h"

#include "HGCSSInfo.hh"
#include "HGCSSEvent.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSParameters.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDigitisation.hh"
#include "HGCSSGenParticle.hh"

#include "HadEnergy.hh"
#include "utilities.h"

#include "TRandom3.h"

using boost::lexical_cast;
namespace po=boost::program_options;

void getValueFromString(std::vector<std::string>& outputString, std::string& sourceString){
     boost::split( outputString, sourceString, boost::is_any_of(","));
};

void getValueFromString(std::vector<int>& outputNumber, std::string& sourceString){
     std::vector<std::string> outputString;
     boost::split( outputString, sourceString, boost::is_any_of(","));
     for(unsigned i(0); i < outputString.size(); i++) outputNumber.push_back(atoi(outputString[i].c_str()));
};

void getValueFromString(std::vector<unsigned>& outputNumber, std::string& sourceString){
     std::vector<std::string> outputString;
     boost::split( outputString, sourceString, boost::is_any_of(","));
     for(unsigned i(0); i < outputString.size(); i++) outputNumber.push_back(atoi(outputString[i].c_str()));
};

void getTotalEnergy(const std::vector<unsigned> & dropLay,
		    TTree* & tree, double & sumEE, double & sumFH, double & sumBH, std::ostringstream& inputsim, std::ostringstream& inputrec, const unsigned nRunsEM=0){

    TFile * simFile = 0;
    TFile * recFile = 0;
 
    TChain  *lSimTree = new TChain("HGCSSTree");
    TChain  *lRecTree = new TChain("RecoTree");

    if (nRunsEM == 0){
      if (!testInputFile(inputsim.str(),simFile)) return;
      lSimTree->AddFile(inputsim.str().c_str());
      if (!testInputFile(inputrec.str(),recFile)) return;
      lRecTree->AddFile(inputrec.str().c_str());
    }
    else {

      for (unsigned i(0);i<nRunsEM;++i){
	std::ostringstream tmpsim,tmprec;	
	tmpsim << inputsim.str() << "_run" << i<< ".root";
	if (!testInputFile(tmpsim.str(),simFile)) continue;
	lSimTree->AddFile(tmpsim.str().c_str());
	tmprec << inputrec.str() << "_run" << i<< ".root";
	if (!testInputFile(tmprec.str(),recFile)) continue;
	lRecTree->AddFile(tmprec.str().c_str());
      }
    }
    bool hasDeDxIn = true;

    HGCSSInfo * info=(HGCSSInfo*)simFile->Get("Info");
    const unsigned versionNumber = info->version();
    std::cout << " -- version number = " << versionNumber << std::endl;
    //if (versionNumber==25 || versionNumber==27 || versionNumber==28) hasDeDxIn=true;

    HGCSSDetector & myDetector = theDetector();
    myDetector.buildDetector(versionNumber,false,false);
    const unsigned nSec = myDetector.nSections();    
    double recSum[nSec];

    for(unsigned iS(0); iS <nSec; iS++){
       recSum[iS] = 0;
    }
 
    HGCSSEvent * event = 0;
    std::vector<HGCSSSamplingSection> * ssvec = 0;
    std::vector<HGCSSSimHit> * simhitvec = 0;
    std::vector<HGCSSRecoHit> * rechitvec = 0;
    std::vector<HGCSSGenParticle> * genvec = 0;

    lSimTree->SetBranchAddress("HGCSSEvent",&event);
    lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
    lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
    lSimTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);

    lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  
    const unsigned nEvts = lSimTree->GetEntries(); 

    std::vector<double> absW;
    bool firstEvent = true;
    for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
      if (ievt%100==0) std::cout << " -- Processing event " << ievt << std::endl;
      lSimTree->GetEntry(ievt);

      if (firstEvent){
	double absweight = 0;
	for (unsigned iL(0); iL<(*ssvec).size();++iL){
	  //double absweight = (*ssvec)[layer].volX0trans()/(*ssvec)[0].volX0trans();
	  bool skipLayer = false;
	  for (unsigned r(0); r<dropLay.size(); r++){
	    if (iL==dropLay[r]) {
	      skipLayer = true;
	      absweight+=hasDeDxIn?(*ssvec)[iL].voldEdx()/(*ssvec)[1].voldEdx() : absWeight(iL,true);
	      break;
	    }
	  }
	  if (skipLayer) {
	    absW.push_back(0);
	    continue;
	  }
	  absweight+=hasDeDxIn?(*ssvec)[iL].voldEdx()/(*ssvec)[1].voldEdx() : absWeight(iL,true);
	  absW.push_back(absweight);
	  absweight = 0;
	}

	std::cout << " -- AbsWeight size: " << absW.size() << std::endl;
      }

      //get truth info
      bool found = false;
      double truthTanx = 0;
      double truthTany = 0;
      //double truthE = 0;
      ROOT::Math::XYZPoint truthPos;
      for (unsigned iP(0); iP<(*genvec).size(); ++iP){//loop on gen particles    
	//if (debug>1 && (*genvec).size()!= 1) 
	//(*genvec)[iP].Print(std::cout);
	if ((*genvec)[iP].trackID()==1){
	  found = true;
	  //set direction
	  truthPos = ROOT::Math::XYZPoint((*genvec)[iP].x(),(*genvec)[iP].y(),(*genvec)[iP].z());
	  truthTanx = (*genvec)[iP].px()/(*genvec)[iP].pz();
	  truthTany = (*genvec)[iP].py()/(*genvec)[iP].pz();
	  //in GeV
	  //truthE = (*genvec)[iP].E()/1000.;
	  //std::cout << " -- true pdgid " << (*genvec)[iP].pdgid()
	  //<< " pos=" << truthPos.x() << " "
	  //<< truthPos.y() << " "
	  //<< truthPos.z() << std::endl;
	  break;
	}
	
      }
      if (!found){
	std::cout << " - Info: no particle G4trackID=1 found, already converted or outside of acceptance ..." << std::endl;
	continue;
      }

      lRecTree->GetEntry(ievt);
      for(unsigned iS(0); iS <nSec; iS++){
        recSum[iS] = 0;
      }
      sumEE = 0;
      sumFH = 0;
      sumBH = 0;

      for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop over rechits
        const HGCSSRecoHit lHit = (*rechitvec)[iH];
    
	double posx = lHit.get_x();
	double posy = lHit.get_y();
	double posz = lHit.get_z();
	double truthposx = truthTanx*(posz-truthPos.z())+truthPos.x();
	double truthposy = truthTany*(posz-truthPos.z())+truthPos.y();
	//consider 5*5 cm^2 section
	if (fabs(posx-truthposx) <= 25 && 
	    fabs(posy-truthposy) <= 25){
	  unsigned layer = lHit.layer();
	  unsigned sec =  myDetector.getSection(layer);
	  double energy = lHit.energy();
	  //std::cout << " -- hit pass: " << energy << " " << absweight << std::endl;
	  
	  recSum[sec] += energy*absW[layer];
	}
      }//loop on hits
    
      for(unsigned iS(0); iS < nSec; iS++){
	//get total energy for each sub-detector
	DetectorEnum type = myDetector.detType(iS); 
	if(type == DetectorEnum::FECAL || type == DetectorEnum::MECAL || type == DetectorEnum::BECAL)
	  sumEE += recSum[iS];
	else if(type == DetectorEnum::FHCAL){
	  sumFH += recSum[iS];
	}
	else if(type == DetectorEnum::BHCAL1 || type == DetectorEnum::BHCAL2)
	  sumBH += recSum[iS];
      }
      //std::cout << " Total E = " << sumEE << " " << sumFH << " " << sumBH << std::endl;
      tree->Fill();
      firstEvent = false;
    }//loop on events
}

void setCalibFactor(const std::vector<unsigned> & dropLay,
		    const std::vector<int> & genEn,
		    const unsigned nRunsEM,
		    double & calibSlopeEE, double & calibOffsetEE,
		    double & calibSlopeFH, double & calibOffsetFH,
		    double & calibSlopeBH, double & calibOffsetBH, 
		    const std::vector<std::pair<std::string, std::string>> & fileName1, 
		    const std::vector<std::pair<std::string, std::string>> & fileName2, 
		    const std::vector<std::pair<std::string, std::string>> & fileName3, 
		    TFile* outputFileEM,
		    TCanvas *c1, TCanvas *c2, TCanvas *c3,
		    TGraphErrors * slopesummaryEE,TGraphErrors * slopesummaryFH,TGraphErrors * slopesummaryBH){

  const unsigned nP = fileName1.size();
  if(fileName1.size() != genEn.size() || fileName1.size() != fileName2.size() || fileName1.size() != fileName3.size()){std::cout << "Got different number of EE/FH/BH files when set calibration factor" << std::endl; return;}
  outputFileEM->cd();
  TTree *treeEE[nP];
  TTree *treeFH[nP];
  TTree *treeBH[nP];

    double E[nP], Energy1[nP], Energy2[nP], Energy3[nP];
    double Energy1err[nP], Energy2err[nP], Energy3err[nP];
    TH1F * p_Energy1=0;
    TH1F * p_Energy2=0;
    TH1F * p_Energy3=0;

    c1->cd();
    c1->Print("EnergyFitsEE.pdf[");
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat("eMRuo");
    //gStyle->SetStatW(0.25);
    //gStyle->SetStatH(0.3);
    std::ostringstream inputsim, inputrec;
    unsigned iP(0);
    for(std::vector<std::pair<std::string, std::string>>::const_iterator iF = fileName1.begin() ; iF != fileName1.end(); ++iF){
      
      inputsim.str("");
      inputrec.str("");
      inputsim << iF->first; 
      inputrec << iF->second;

      E[iP] = genEn[iP];
      std::ostringstream var;
      var << "treeEE_" << E[iP];
      outputFileEM->cd();
      treeEE[iP] = new TTree(var.str().c_str(),"EE energies");
      double sumEE(0), sumFH(0), sumBH(0);   
      var.str("");
      var << "E_EE";
      treeEE[iP]->Branch(var.str().c_str(),&sumEE);
      var.str("");
      var << "E_FH";
      treeEE[iP]->Branch(var.str().c_str(),&sumFH);
      var.str("");
      var << "E_BH";
      treeEE[iP]->Branch(var.str().c_str(),&sumBH);

      getTotalEnergy(dropLay,treeEE[iP], sumEE, sumFH, sumBH, inputsim, inputrec,nRunsEM);
      //for(unsigned iE(0); iE < Etotal[0].size(); iE++)p_Energy1->Fill(Etotal[0][iE]);
      outputFileEM->cd();
      treeEE[iP]->Write();
      //std::cout << " -- Entries in tree after fill: " << treeEE[iP]->GetEntries() << std::endl;

      //treeEE[iP]->Show(0);
      var.str("");
      var << "E_EE";
      treeEE[iP]->Draw(var.str().c_str());
      p_Energy1 = (TH1F*)(gPad->GetPrimitive("htemp"))->Clone("p_Energy1");
      p_Energy1->Draw("PE");
      TF1 *fitResult = new TF1("fitResult","[0]*TMath::Gaus(x,[1],[2],0)",p_Energy1->GetXaxis()->GetXmin(),p_Energy1->GetXaxis()->GetXmax());
      fitResult->SetParameters(p_Energy1->Integral(),
			       p_Energy1->GetXaxis()->GetBinCenter(p_Energy1->GetMaximumBin()),
			       fabs(p_Energy1->GetRMS()));
      p_Energy1->Fit("fitResult","L0QEMI","",
		     fitResult->GetParameter(1)-2*fitResult->GetParameter(2),
		     fitResult->GetParameter(1)+2*fitResult->GetParameter(2));
      fitResult->SetLineColor(2);
      fitResult->Draw("same");
      //c1->Print("checkE3.png");
      c1->Print("EnergyFitsEE.pdf");

      Energy1[iP] = fitResult->GetParameter(1);//p_Energy1->GetMean();
      Energy1err[iP] = fitResult->GetParError(1);//p_Energy1->GetMeanError();

      std::cout << "Total energy  gen=" << E[iP] << " rec= " << Energy1[iP] << "(MIPS)" << std::endl;
      p_Energy1->Delete();
      iP++;
    } 
    c1->Print("EnergyFitsEE.pdf]");

    c2->cd();
    c2->Print("EnergyFitsFH.pdf[");

    iP = 0;
    for(std::vector<std::pair<std::string, std::string>>::const_iterator iF = fileName2.begin() ; iF != fileName2.end(); ++iF){
      
      inputsim.str("");
      inputrec.str("");
      inputsim << iF->first; 
      inputrec << iF->second;
     
      std::ostringstream var;
      var << "treeFH_" << E[iP];
      outputFileEM->cd();
      treeFH[iP] = new TTree(var.str().c_str(),"FH+BH energies");
      double sumEE(0), sumFH(0), sumBH(0);   
      var.str("");
      var << "E_EE";
      treeFH[iP]->Branch(var.str().c_str(),&sumEE);
      var.str("");
      var << "E_FH";
      treeFH[iP]->Branch(var.str().c_str(),&sumFH);
      var.str("");
      var << "E_BH";
      treeFH[iP]->Branch(var.str().c_str(),&sumBH);
      getTotalEnergy(dropLay,treeFH[iP], sumEE, sumFH, sumBH, inputsim, inputrec,nRunsEM);
      //for(unsigned iE(0); iE < Etotal[1].size(); iE++)p_Energy2->Fill(Etotal[1][iE]);
      outputFileEM->cd();
      treeFH[iP]->Write();
      var.str("");
      var << "E_FH";
      treeFH[iP]->Draw(var.str().c_str());
      p_Energy2 = (TH1F*)(gPad->GetPrimitive("htemp"))->Clone("p_Energy2");
      p_Energy2->Draw("PE");
      TF1 *fitResult = new TF1("fitResult","[0]*TMath::Gaus(x,[1],[2],0)",p_Energy2->GetXaxis()->GetXmin(),p_Energy2->GetXaxis()->GetXmax());
      fitResult->SetParameters(p_Energy2->Integral(),
			       p_Energy2->GetXaxis()->GetBinCenter(p_Energy2->GetMaximumBin()),
			       p_Energy2->GetRMS());
      p_Energy2->Fit("fitResult","L0QEMI","",
		     fitResult->GetParameter(1)-2*fitResult->GetParameter(2),
		     fitResult->GetParameter(1)+2*fitResult->GetParameter(2));
      fitResult->SetLineColor(2);
      fitResult->Draw("same");

      c2->Print("EnergyFitsFH.pdf");

      Energy2[iP] = fitResult->GetParameter(1);//p_Energy2->GetMean();
      Energy2err[iP] = fitResult->GetParError(1);//p_Energy2->GetMeanError();


      std::cout << "Total energy  gen=" << E[iP] << " rec  = " << Energy2[iP] << "(MIPS)" << std::endl;
      p_Energy2->Delete();
      iP++;
    } 
    c2->Print("EnergyFitsFH.pdf]");
    
    c3->cd();
    c3->Print("EnergyFitsBH.pdf[");

    iP = 0;
    for(std::vector<std::pair<std::string, std::string>>::const_iterator iF = fileName3.begin() ; iF != fileName3.end(); ++iF){
      
      inputsim.str("");
      inputrec.str("");
      inputsim << iF->first; 
      inputrec << iF->second;
     
      std::ostringstream var;
      var << "treeBH_" << E[iP];
      outputFileEM->cd();
      treeBH[iP] = new TTree(var.str().c_str(),"BH energies");
      double sumEE(0), sumFH(0), sumBH(0);   
      var.str("");
      var << "E_EE";
      treeBH[iP]->Branch(var.str().c_str(),&sumEE);
      var.str("");
      var << "E_FH";
      treeBH[iP]->Branch(var.str().c_str(),&sumFH);
      var.str("");
      var << "E_BH";
      treeBH[iP]->Branch(var.str().c_str(),&sumBH);
      getTotalEnergy(dropLay,treeBH[iP], sumEE, sumFH, sumBH, inputsim, inputrec,nRunsEM);
      outputFileEM->cd();
      treeBH[iP]->Write();
     //for(unsigned iE(0); iE < Etotal[1].size(); iE++)p_Energy3->Fill(Etotal[2][iE]);
      var.str("");
      var << "E_BH";
      treeBH[iP]->Draw(var.str().c_str());
      p_Energy3 = (TH1F*)(gPad->GetPrimitive("htemp"))->Clone("p_Energy3");

      p_Energy3->Draw("PE");
      TF1 *fitResult = new TF1("fitResult","[0]*TMath::Gaus(x,[1],[2],0)",p_Energy3->GetXaxis()->GetXmin(),p_Energy3->GetXaxis()->GetXmax());
      fitResult->SetParameters(p_Energy3->Integral(),
			       p_Energy3->GetXaxis()->GetBinCenter(p_Energy3->GetMaximumBin()),
			       p_Energy3->GetRMS());
      p_Energy3->Fit("fitResult","L0QEMI","",
		     fitResult->GetParameter(1)-2*fitResult->GetParameter(2),
		     fitResult->GetParameter(1)+2*fitResult->GetParameter(2));
      fitResult->SetLineColor(2);
      fitResult->Draw("same");

      c3->Print("EnergyFitsBH.pdf");

      Energy3[iP] = fitResult->GetParameter(1);//p_Energy3->GetMean();
      Energy3err[iP] = fitResult->GetParError(1);//p_Energy3->GetMeanError();

      std::cout << "Total energy  gen=" << E[iP] << " rec  = " << Energy3[iP] << "(MIPS)" << std::endl;
      p_Energy3->Delete();
      iP++;
    }
    c3->Print("EnergyFitsBH.pdf]");

    for(unsigned iP(0); iP < nP; iP++){
      slopesummaryEE->SetPoint(iP,E[iP],Energy1[iP]);
      slopesummaryEE->SetPointError(iP,0,Energy1err[iP]);
      slopesummaryFH->SetPoint(iP,E[iP],Energy2[iP]);
      slopesummaryFH->SetPointError(iP,0,Energy2err[iP]);
      slopesummaryBH->SetPoint(iP,E[iP],Energy3[iP]);
      slopesummaryBH->SetPointError(iP,0,Energy3err[iP]);
    }
    /*    c1->cd();
    slopesummaryEE->SetMarkerSize(1);
    slopesummaryEE->SetMarkerStyle(20);
    slopesummaryEE->SetMarkerColor(1);
    slopesummaryEE->Draw("AP");
    slopesummaryEE->Fit("pol1");
    TF1 *fitEE = (TF1*)slopesummaryEE->GetFunction("pol1");
   // fit->Draw("same");
    calibOffsetEE = fitEE->GetParameter(0);
    calibSlopeEE = fitEE->GetParameter(1);

    c2->cd();
    slopesummaryFH->SetMarkerSize(1);
    slopesummaryFH->SetMarkerStyle(20);
    slopesummaryFH->SetMarkerColor(1);
    slopesummaryFH->Draw("AP");
    slopesummaryFH->Fit("pol1");
    TF1 *fitFH = (TF1*)slopesummaryFH->GetFunction("pol1");
   // fit->Draw("same");
    calibOffsetFH = fitFH->GetParameter(0);
    calibSlopeFH = fitFH->GetParameter(1);

    c3->cd();
    slopesummaryBH->SetMarkerSize(1);
    slopesummaryBH->SetMarkerStyle(20);
    slopesummaryBH->SetMarkerColor(1);
    slopesummaryBH->Draw("AP");
    slopesummaryBH->Fit("pol1");
    TF1 *fitBH = (TF1*)slopesummaryBH->GetFunction("pol1");
   // fit->Draw("same");
    calibOffsetBH = fitBH->GetParameter(0);
    calibSlopeBH = fitBH->GetParameter(1);
    */

    return;
}

/*void setCalibFactor(double & calibFactor, const std::vector<std::pair<std::string, std::string>>& fileName,const double & FHtoEslope,const double & FHtoEoffset,const double & BHtoEslope,const double & BHtoEoffset,  TCanvas *c){

    c->cd();
    TGraph *slopesummary = new TGraph();
    //const unsigned nP = fileName.size();
    TH2F *p_FHvsBH = new TH2F("","",1000,0,10000,1000,0,10000);
    std::vector<std::vector<double>> Etotal;

    std::ostringstream inputsim, inputrec;
    for(std::vector<std::pair<std::string, std::string>>::const_iterator iF = fileName.begin() ; iF != fileName.end(); ++iF){
      
      inputsim.str("");
      inputrec.str("");
      inputsim << iF->first; 
      inputrec << iF->second;
     
      getTotalEnergy(Etotal, inputsim, inputrec);
      unsigned nevt = Etotal[1].size();
      for(unsigned ievt(0); ievt < nevt; ievt++){
         p_FHvsBH->Fill(Etotal[0][ievt]+(Etotal[1][ievt]-FHtoEoffset)/FHtoEslope, (Etotal[2][ievt]-BHtoEoffset)/BHtoEslope);
      }

      TProfile *prof = p_FHvsBH->ProfileX();
      prof->Fit("pol1");
      TF1 *f1 = (TF1*)prof->GetFunction("pol1");
      double slope = f1->GetParameter(1);
  
      Int_t np=slopesummary->GetN();  
      slopesummary->SetPoint(np,np,slope);
    }
 
    slopesummary->Draw("AP"); 
    double slopeMean = fabs(slopesummary->GetMean(2));
    calibFactor = slopeMean;
    } */
         

int main(int argc, char** argv){//main  

  //Input output and config options
  std::string cfg;
  bool concept;
  //size of signal region to perform Chi2 position fit.
  unsigned pNevts;
  bool doFHEMcalib;
  bool doBHEMcalib;
  bool doFHBHcalib;
  std::string eeSimFiles;
  std::string fhSimFiles;
  std::string bhSimFiles;
  std::string eeRecoFiles;
  std::string fhRecoFiles;
  std::string bhRecoFiles;
  std::string filePath;
  std::string digifilePath;
  unsigned nRuns;
  unsigned nRunsEM;
  std::string genEnergy;
  std::string genEnergyEM;
  std::string simFileName;
  std::string recoFileName;
  std::string outPath;
  std::string outFileEM;
  std::string outFilePi;
  unsigned nSiLayers;
  //0:do just the energies, 1:do fit+energies, 2: do zpos+fit+energies
  unsigned debug;
  std::string dropLayers;

  po::options_description preconfig("Configuration"); 
  preconfig.add_options()("cfg,c",po::value<std::string>(&cfg)->required());
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(preconfig).allow_unregistered().run(), vm);
  po::notify(vm);
  po::options_description config("Configuration");
  config.add_options()
    //Input output and config options //->required()
    ("concept",        po::value<bool>(&concept)->default_value(true))
    ("pNevts,n",       po::value<unsigned>(&pNevts)->default_value(0))
    ("doFHEMcalib", po::value<bool>(&doFHEMcalib)->default_value(false))
    ("doBHEMcalib", po::value<bool>(&doBHEMcalib)->default_value(false))
    ("doFHBHcalib",  po::value<bool>(&doFHBHcalib)->default_value(false))
    ("eeSimFiles",     po::value<std::string>(&eeSimFiles)->default_value(""))
    ("fhSimFiles",     po::value<std::string>(&fhSimFiles)->default_value(""))
    ("bhSimFiles",     po::value<std::string>(&bhSimFiles)->default_value(""))
    ("eeRecoFiles",     po::value<std::string>(&eeRecoFiles)->default_value(""))
    ("fhRecoFiles",     po::value<std::string>(&fhRecoFiles)->default_value(""))
    ("bhRecoFiles",     po::value<std::string>(&bhRecoFiles)->default_value(""))
    ("filePath,i",     po::value<std::string>(&filePath)->required())
    ("nRuns",        po::value<unsigned>(&nRuns)->default_value(0))
    ("nRunsEM",        po::value<unsigned>(&nRunsEM)->default_value(0))
    ("genEnergy",    po::value<std::string>(&genEnergy)->required()) 
    ("genEnergyEM",    po::value<std::string>(&genEnergyEM)->required()) 
    ("simFileName,s",  po::value<std::string>(&simFileName)->required())
    ("recoFileName,r", po::value<std::string>(&recoFileName)->required())
    ("outPath,o",      po::value<std::string>(&outPath)->required())
    ("outFileEM,o",      po::value<std::string>(&outFileEM)->required())
    ("outFilePi,o",      po::value<std::string>(&outFilePi)->required())
    ("nSiLayers",      po::value<unsigned>(&nSiLayers)->default_value(2))
    ("debug,d",        po::value<unsigned>(&debug)->default_value(0))
    ("dropLayers",    po::value<std::string>(&dropLayers)->required()) 
    ;

  // ("output_name,o",            po::value<std::string>(&outputname)->default_value("tmp.root"))

  po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
  po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
  po::notify(vm);

  //////////////////////////////////////////////////////////
  //// Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////
//  bool isCalibed = true;
  //bool isCalibed = false;
  //////////////////////////////////////////////////////////
  //// Hardcoded factor ////////////////////////////////////
  //////////////////////////////////////////////////////////
  double FHtoEslope = 21.02;
  double FHtoEoffset = 4;
  double BHtoEslope = 8.62;
  double BHtoEoffset = 3.7;
  double ECALslope = 92.83;
  double ECALoffset = 13;

  double FHtoBHslope = 1.49;
  double EEtoHslope = 5.28;
  double EEtoHslopePar0 = 30.;
  double EEtoHslopePar1 = 5.31;


  //double HcalPionOffset = 0.; 
  //double HcalPionCalib = 0.92;
////////////////////////////////////
//


  std::vector<unsigned> dropLay;
  getValueFromString(dropLay, dropLayers);

  const unsigned nDrop = dropLay.size();
  std::cout << " -- Dropping " << nDrop << " layers: " << dropLayers << std::endl;

  if(doFHEMcalib && doBHEMcalib){

    if(eeSimFiles == ""|| eeRecoFiles == "" || fhSimFiles == "" || fhRecoFiles == ""|| bhSimFiles == "" || bhRecoFiles == "") {std::cout << "EE/FH/BH files cannot be empty for FHvsEE and BHvsEE calibration" << std::endl;}
    else{
      std::vector<std::string> eeSimFileVec, eeRecoFileVec, fhSimFileVec, fhRecoFileVec, bhSimFileVec, bhRecoFileVec;
      getValueFromString( eeSimFileVec, eeSimFiles);
      getValueFromString( eeRecoFileVec,eeRecoFiles);
      getValueFromString( fhSimFileVec,fhSimFiles);
      getValueFromString( fhRecoFileVec,fhRecoFiles);
      getValueFromString( bhSimFileVec,bhSimFiles);
      getValueFromString( bhRecoFileVec,bhRecoFiles);
      std::vector<std::pair<std::string, std::string>> eeFileName, fhFileName, bhFileName;
      for(unsigned iF(0); iF<eeSimFileVec.size(); iF++){
         eeFileName.push_back(make_pair(eeSimFileVec[iF], eeRecoFileVec[iF]));
         fhFileName.push_back(make_pair(fhSimFileVec[iF], fhRecoFileVec[iF]));
         bhFileName.push_back(make_pair(bhSimFileVec[iF], bhRecoFileVec[iF]));
      }
      TCanvas *c1 = new TCanvas("EEcalibtoEM","EEcalibtoEM",800,600);
      TCanvas *c2 = new TCanvas("FHcalibtoEM","FHcalibtoEM",800,600);
      TCanvas *c3 = new TCanvas("BHcalibtoEM","BHcalibtoEM",800,600);

      std::vector<int> genEn;
      getValueFromString(genEn, genEnergyEM);

    /////////////////////////////////////////////////////////////
    //output
    ///////////////////////////////////////////////////////////////
    std::ostringstream outFile;
    outFile.str("");
    if (outFileEM.find(".root") == outFileEM.npos)
      outFile << outPath << "/" << outFileEM << "_" << genEn[0] << ".root";
    else 
      outFile << outPath << "/" << outFileEM;
    TFile *outputFileEM = TFile::Open(outFile.str().c_str(),"RECREATE");

    if (!outputFileEM) {
      std::cout << " -- Error, output file " << outFile.str() << " cannot be opened. Please create output directory. Exiting..." << std::endl;
      return 1;
    }
     else {
       std::cout << " -- output file " << outputFileEM->GetName() << " successfully opened." << std::endl;
     }
    outputFileEM->cd();
    TGraphErrors * slopesummaryEE = new TGraphErrors();
    TGraphErrors * slopesummaryFH = new TGraphErrors();
    TGraphErrors * slopesummaryBH = new TGraphErrors();


    setCalibFactor(dropLay,genEn,nRunsEM,
		     ECALslope,ECALoffset,
		     FHtoEslope,FHtoEoffset,
		     BHtoEslope,BHtoEoffset,
		     eeFileName, fhFileName, bhFileName,
		     outputFileEM,
		     c1,c2,c3,
		     slopesummaryEE,slopesummaryFH,slopesummaryBH); 
      //c1->SaveAs("EEcalibtoEM.png");
      //c2->SaveAs("FHcalibtoEM.png");
      //c3->SaveAs("BHcalibtoEM.png");

      //std::cout << " Summary of calibration constant to EM scale: " << std::endl
      //<< " ECAL : slope = " << ECALslope << " offset = " << ECALoffset << std::endl
      //<< " FHCAL: slope = " << FHtoEslope << " offset = " << FHtoEoffset << std::endl
      //<< " BHCAL: slope = " << BHtoEslope << " offset = " << BHtoEoffset << std::endl;
      //outputFileEM->cd();
      //slopesummaryEE->Write();
      //slopesummaryFH->Write();
      //slopesummaryBH->Write();
      outputFileEM->Close();
    }
    return 0;
  }
  /*  if(doBHEMcalib){
    if(eeSimFiles == ""|| eeRecoFiles == "" || bhSimFiles == "" || bhRecoFiles == "") {std::cout << "EE/FH files cannot be empty for FHvsEE calibration" << std::endl;}
    else{
      std::vector<std::string> eeSimFileVec, eeRecoFileVec, bhSimFileVec, bhRecoFileVec;
      getValueFromString( eeSimFileVec,eeSimFiles);
      getValueFromString( eeRecoFileVec,eeRecoFiles);
      getValueFromString( bhSimFileVec,bhSimFiles);
      getValueFromString( bhRecoFileVec,bhRecoFiles);
      std::vector<std::pair<std::string, std::string>> eeFileName, bhFileName;
      for(unsigned iF(0); iF<eeSimFileVec.size(); iF++){
         eeFileName.push_back(make_pair(eeSimFileVec[iF], eeRecoFileVec[iF]));
         bhFileName.push_back(make_pair(bhSimFileVec[iF], bhRecoFileVec[iF]));
      }
      TCanvas *c2 = new TCanvas("BHcalibtoEM","BHcalibtoEM",800,600);
      setCalibFactor(BHtoEslope,BHtoEoffset,eeFileName, bhFileName, c2); 
      c2->SaveAs("BHcalibtoEM.png");
    }
    }
//this is to be done on pion files...
  if(doFHBHcalib){
    if(fhSimFiles == ""|| fhRecoFiles == "") {std::cout << "FH files cannot be empty for FHvsEE calibration" << std::endl;}
    else{
      std::vector<std::string> fhSimFileVec, fhRecoFileVec;
      getValueFromString( fhSimFileVec,fhSimFiles);
      getValueFromString( fhRecoFileVec,fhRecoFiles);
      std::vector<std::pair<std::string, std::string>> fhFileName;
      for(unsigned iF(0); iF<fhSimFileVec.size(); iF++){
         fhFileName.push_back(make_pair(fhSimFileVec[iF], fhRecoFileVec[iF]));
      }
      TCanvas *c3 = new TCanvas("BHcalibtoFH","BHcalibtoFH",800,600);
      setCalibFactor(FHtoBHslope, fhFileName, FHtoEslope,FHtoEoffset,BHtoEslope,BHtoEoffset,c3); 
    }
    }*/
      

  //////////////////////////////////////////////////////////
  //// End Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////

  TRandom3 *lRndm = new TRandom3();
  lRndm->SetSeed(1234);

  std::cout << " -- Random3 seed = " << lRndm->GetSeed() << std::endl
	    << " ----------------------------------------" << std::endl;

  std::string inFilePath = filePath+simFileName;

  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << filePath << std::endl
	    << " -- Output file path: " << outPath << std::endl
	    << " -- Requiring " << nSiLayers << " si layers." << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events." << std::endl;

  std::vector<int> genEn;
  getValueFromString(genEn, genEnergy);
  const unsigned nGenEn=genEn.size();
  //double eta = 2.00;

  /////////////////////////////////////////////////////////////
  //input file format
  ///////////////////////////////////////////////////////////////

  std::size_t begin = simFileName.find("_et")+3;
  std::size_t middle = simFileName.find("_eta");
  std::size_t end = simFileName.find_last_of(".root")+1;
  std::size_t run(0);
  if(nRuns>0)run = simFileName.find("_run");
  std::string simHeader = simFileName.substr(0,begin);
  std::string simAppend("");
  if(nRuns>0)simAppend = simFileName.substr(middle,run-middle);
  else simAppend = simFileName.substr(middle,end-middle);

  std::size_t begin_rec = recoFileName.find("_et")+3;
  std::size_t middle_rec = recoFileName.find("_eta");
  std::size_t end_rec = recoFileName.find_last_of(".root")+1;
  std::size_t run_rec(0);
  if(nRuns>0)run_rec = recoFileName.find("_run");
  std::string recHeader = recoFileName.substr(0,begin_rec);
  std::string recAppend;
  if(nRuns>0)recAppend = recoFileName.substr(middle_rec,run_rec-middle_rec);
  else recAppend = recoFileName.substr(middle_rec,end_rec-middle_rec);


  //std::vector<double> EE[nGenEn];
  //std::vector<double> EFHCAL[nGenEn];
  //std::vector<double> EBHCAL[nGenEn];
  //std::vector<double> GlobalC[nGenEn];

  for(unsigned iGen(0); iGen < nGenEn; iGen++){
 
     TFile * simFile = 0;
     TFile * recFile = 0;
  
     TChain  *lSimTree = new TChain("HGCSSTree");
     TChain  *lRecTree = new TChain("RecoTree");
   
     std::ostringstream inputsim, inputrec;
     if (nRuns == 0){
       inputsim << filePath << "/" << simHeader << genEn[iGen] << simAppend;
       inputrec << filePath << "/" << recHeader << genEn[iGen] << recAppend;;
       if (!testInputFile(inputsim.str(),simFile)) return 1;
       lSimTree->AddFile(inputsim.str().c_str());
       if (!testInputFile(inputrec.str(),recFile)) return 1;
       lRecTree->AddFile(inputrec.str().c_str());
     }
     else {
       for (unsigned i(0);i<nRuns;++i){
         inputsim.str("");
         inputsim << filePath << "/" << simHeader << genEn[iGen] << simAppend << "_run" << i << ".root";
         inputrec.str("");
         inputrec << filePath << "/" << recHeader << genEn[iGen] << recAppend << "_run" << i << ".root";
         if (!testInputFile(inputsim.str(),simFile)) continue;
         lSimTree->AddFile(inputsim.str().c_str());
         if (!testInputFile(inputrec.str(),recFile)) continue;
         lRecTree->AddFile(inputrec.str().c_str());
       }
     }


     /////////////////////////////////////////////////////////////
     //output
     ///////////////////////////////////////////////////////////////
     std::ostringstream outFile;
     outFile.str("");
     outFile << outPath << outFilePi << "_et" << genEn[iGen] << ".root";
     TFile *outputFile = TFile::Open(outFile.str().c_str(),"RECREATE");

     if (!outputFile) {
       std::cout << " -- Error, output file " << outFile.str() << " cannot be opened. Please create output directory. Exiting..." << std::endl;
     return 1;
     }
     else {
       std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
     }
   
     ///////////////////////////////
     // Info    ////////////////////
     ///////////////////////////////
   
     HGCSSInfo * info=(HGCSSInfo*)simFile->Get("Info");
     const double cellSize = info->cellSize();
     const unsigned versionNumber = info->version();
     const unsigned model = info->model();

     bool isCaliceHcal = versionNumber==23;   
     std::cout << " -- Version number is : " << versionNumber 
   	    << ", model = " << model
   	    << ", cellSize = " << cellSize
   	    << std::endl;

     //initialise detector
     HGCSSDetector & myDetector = theDetector();
 
     myDetector.buildDetector(versionNumber,concept,isCaliceHcal);

     const unsigned nLayers = myDetector.nLayers();
     const unsigned nSections = myDetector.nSections();

     std::cout << " -- N layers = " << nLayers << std::endl
      	       << " -- N sections = " << nSections << std::endl;
     
     HadEnergy myhadReso(dropLay,myDetector, lSimTree, lRecTree, outputFile, pNevts);
     myhadReso.addLimMIP(3);
     myhadReso.addLimMIP(5);
     myhadReso.addLimMIP(10);
     myhadReso.addLimMIP(15);
     myhadReso.addLimMIP(20);
     myhadReso.addLimMIP(25);
     myhadReso.addLimMIP(30);
     myhadReso.addLimMIP(35);
     myhadReso.addLimMIP(40);
     myhadReso.addLimMIP(50);
     myhadReso.bookHist(outputFile);
     myhadReso.setFHtoE(FHtoEslope, FHtoEoffset);
     myhadReso.setBHtoE(BHtoEslope, BHtoEoffset);
     myhadReso.setEEcalib(ECALslope, ECALoffset);
     myhadReso.setFHtoBH(FHtoBHslope);
     myhadReso.setEEtoH(EEtoHslope);
     myhadReso.setEEtoHPar(EEtoHslopePar0,EEtoHslopePar1);
 
     myhadReso.fillEnergies(); 

     outputFile->Close();
  }// loop over GenEn 
  
  return 0;
  
  
}//main
