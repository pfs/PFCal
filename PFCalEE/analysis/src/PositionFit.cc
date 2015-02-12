#include <iomanip>

#include "PositionFit.hh"
#include "Clusterizer.hh"
#include "HGCSSEvent.hh"
#include "HGCSSInfo.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSDetector.hh"

double PositionFit::getW0(const unsigned layer){
  if (layer==0) return 1.6;
  if (layer==1) return 2.15;
  if (layer==2) return 2.6;
  if (layer==3) return 2.4;
  if (layer==4) return 2.5;
  if (layer==5) return 3.15;
  if (layer==6) return 3.3;
  if (layer==7) return 2.9;
  if (layer==8) return 3.1;
  if (layer==9) return 2.6;
  if (layer==10) return 2.9;
  if (layer==11) return 2.5;
  if (layer==12) return 2.7;
  if (layer==13) return 2.3;
  if (layer==14) return 2.5;
  if (layer==15) return 2.1;
  if (layer==16) return 2.3;
  if (layer==17) return 1.95;
  if (layer==18) return 2.15;
  if (layer==19) return 1.9;
  if (layer==20) return 2.05;
  if (layer==21) return 1.8;
  if (layer==22) return 2.15;
  if (layer==23) return 1.9;
  if (layer==24) return 2.1;
  if (layer==25) return 1.8;
  if (layer==26) return 1.9;
  if (layer==27) return 3;
  if (layer==28) return 3;
  if (layer==29) return 3;

  return 0;
}

double PositionFit::DeltaPhi(const double & phi1, const double & phi2){
  double dphi = phi1 - phi2;
  if (dphi< (-1.*TMath::Pi())) dphi += 2*TMath::Pi();
  if (dphi>TMath::Pi()) dphi -= 2*TMath::Pi();
  return dphi;
}

std::pair<unsigned, std::pair<double,double> > PositionFit::findMajorityValue(std::vector<std::pair<double,double> > & values) const
{
  
  unsigned lTot = values.size();
  if (!lTot) return std::pair<unsigned, std::pair<double,double> >(0,std::pair<double,double>(0,0));
  
  std::sort(values.begin(),values.end());

  unsigned lMajorityCounter = 0;
  std::pair<double,double> lMaj = std::pair<double,double>(0,0);
  
  std::vector<std::pair<double,double> >::iterator lIter = values.begin();
  for ( ; lIter != values.end(); ) {
    //std::cout << lIter->first << " " << lIter->second << std::endl;
    unsigned lCounter = std::count(lIter,values.end(),*lIter);
    if (lCounter > lMajorityCounter) {
      lMajorityCounter = lCounter;
      lMaj = *lIter;
    }
    lIter += lCounter;
  }
    
  //std::cout << " -- Found majority value " << lMaj.first << " " << lMaj.second << " for " << lMajorityCounter << " elements out of " << values.size() << "." << std::endl;

  return std::pair<unsigned, std::pair<double,double> >(lMajorityCounter,lMaj);
  
}

PositionFit::PositionFit(const unsigned nSR,
			 const double & residualMax, 
			 const unsigned nLayers, 
			 const unsigned nSiLayers,
			 const bool applyPuMixFix,
			 const unsigned debug,
			 const bool doMatrix
			 ){
  nSR_ = nSR;
  residualMax_ = residualMax;
  chi2ndfmax_ = 10;
  precision_ = 1.e-10;
  seedMipThreshold_ = 10;
  maxdR_ = 0.1;
  nLayers_ = nLayers;
  nSiLayers_ = nSiLayers;
  debug_ = debug;
  useMeanPU_ = true;
  fixForPuMixBug_ = applyPuMixFix;
  doMatrix_ = doMatrix;
  saveEtree_ = true;
  doLogWeight_ = true;

  p_nGenParticles = 0;
  //p_numberOfMaxTried = 0;
  //p_dRMaxTruth = 0;
  p_hitMeanPuContrib = 0;
  p_hitEventPuContrib = 0;
  p_diffPuContrib = 0;

  p_etavsphi = 0;
  p_etavsphi_max = 0;
  p_etavsphi_truth = 0;
  p_yvsx_max = 0;
  p_yvsx_truth = 0;

  p_residuals_x = 0;
  p_residuals_y = 0;
  p_errorMatrix_xy = 0;
  p_corrMatrix_xy = 0;
  //p_errorMatrix_y = 0;
  //p_corrMatrix_y = 0;
  p_chi2[0] = 0;
  p_chi2[1] = 0;
  p_iterations = 0;

  p_chi2overNDF[0] = 0;
  p_impactX0[0] = 0;
  p_impactY0[0] = 0;
  p_impactX0_residual = 0;
  p_impactY0_residual = 0;
  p_impactXFF[0] = 0;
  p_impactYFF[0] = 0;
  p_impactXFF_residual = 0;
  p_impactYFF_residual = 0;
  p_impactX14[0] = 0;
  p_impactY14[0] = 0;
  p_impactX14_residual = 0;
  p_impactY14_residual = 0;
  p_impactZx[0] = 0;
  p_impactZx[1] = 0;
  p_impactZx_residual = 0;
  p_impactZy[0] = 0;
  p_impactZy[1] = 0;
  p_impactZy_residual = 0;

  p_tanAngleX[0] = 0;
  p_tanAngleY[0] = 0;
  p_tanAngleX_residual = 0;
  p_tanAngleY_residual = 0;
  p_angleX_residual = 0;
  p_angleY_residual = 0;
  p_eta_reco = 0;
  p_phi_reco = 0;
  p_eta_truth = 0;
  p_phi_truth = 0;
  p_eta_residual = 0;
  p_phi_residual = 0;
  p_positionReso = 0;
  p_angularReso = 0;
  p_chi2overNDF[1] = 0;
  p_impactX0[1] = 0;
  p_impactY0[1] = 0;
  p_impactXFF[1] = 0;
  p_impactYFF[1] = 0;
  p_impactX14[1] = 0;
  p_impactY14[1] = 0;
  p_tanAngleX[1] = 0;
  p_tanAngleY[1] = 0;

}

void PositionFit::initialise(TFile *outputFile,
			     const std::string outputDir,
			     const std::string outFolder, 
			     const HGCSSGeometryConversion & geomConv, 
			     const HGCSSPUenergy & puDensity){

  outputDir_ = outputDir;
  setOutputFile(outputFile);

  outFolder_ = outFolder;
  matrixFolder_ = outFolder_;

  geomConv_ = geomConv;
  puDensity_ = puDensity;

  nL_mean_.resize(nLayers_*2,0);
  mean_.resize(nLayers_*2,0);

  nL_sigma_.resize(nLayers_*2,nL_mean_);
  sigma_.resize(nLayers_*2,mean_);
  for (unsigned iL(0);iL<2*nLayers_;++iL){
    nL_sigma_[iL].resize(nLayers_*2,0);
    sigma_[iL].resize(nLayers_*2,0);
  }

}

void PositionFit::initialiseClusterHistograms(){

  outputFile_->cd(outputDir_.c_str());

  p_nClusters = new TH1F("p_nClusters",";n_{clusters};n_{events}",1000,0,3000);

  p_clusnHits_all = new TH1F("p_clusnHits_all",";n_{hits};n_{clusters}",500,0,1000);
  p_seedEoverE_all = new TH1F("p_seedEoverE_all",";seedE/E;n_{clusters}",100,0,1);
  p_clusLayer_all = new TH1F("p_clusLayer_all",";cluster layer;n_{clusters}",nLayers_,0,nLayers_);
  p_clusWidth_all = new TH1F("p_clusWidth_all",";cluster width (layers);n_{clusters}",nLayers_,0,nLayers_);
  p_seeddeta_all = new TH1F("p_seeddeta_all",";#Delta#eta(seed,cluster);n_{clusters}",100,-0.5,0.5);
  p_seeddphi_all = new TH1F("p_seeddphi_all",";#Delta#phi(seed,cluster);n_{clusters}",100,-0.5,0.5);

  p_clusnHits_sel = new TH1F("p_clusnHits_sel",";n_{hits};n_{events}",500,0,1000);
  p_seedEoverE_sel = new TH1F("p_seedEoverE_sel",";seedE/E;n_{events}",100,0,1);
  p_clusLayer_sel = new TH1F("p_clusLayer_sel",";cluster layer;n_{events}",nLayers_,0,nLayers_);
  p_clusWidth_sel = new TH1F("p_clusWidth_sel",";cluster width (layers);n_{events}",nLayers_,0,nLayers_);
  p_seeddeta_sel = new TH1F("p_seeddeta_sel",";#Delta#eta(seed,cluster);n_{events}",100,-0.1,0.1);
  p_seeddphi_sel = new TH1F("p_seeddphi_sel",";#Delta#phi(seed,cluster);n_{events}",100,-0.1,0.1);

   p_mindRtruth = new TH1F("p_mindRtruth",";mindR(cluster,truth);n_{events}",100,0,0.5);

 }

 void PositionFit::initialisePositionHistograms(){
   //for yvsx plots
   unsigned nX=2*1700/10-2;
   double minX=-1.*nX*5-5,maxX=nX*5+5;
   double minY=minX,maxY=maxX;
   nX += 1;
   unsigned nY = nX;

   outputFile_->cd(outputDir_.c_str());

   p_nGenParticles = new TH1F("p_nGenParticles",";nGenParticles",10,0,10);
   //p_numberOfMaxTried = new TH1F("p_numberOfMaxTried",";max cells tried",250,0,250);
   //p_dRMaxTruth = new TH1F("p_dRMaxTruth",";#Delta R(max,truth)",200,0,0.1);
   //p_dRMaxTruth->StatOverflows();

   p_genxy.resize(nLayers_,0);
   p_recoxy.resize(nLayers_,0);
   p_dRmin.resize(nLayers_,0);
   std::ostringstream lName;
   for (unsigned iL(0); iL<nLayers_; ++iL){
     lName.str("");
     lName << "p_genxy_" << iL;
     p_genxy[iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",
			    nX*10,minX,maxX,
			    nY*10,minY,maxY);
     lName.str("");
     lName << "p_recoxy_" << iL;
     p_recoxy[iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",
			     nX,minX,maxX,
			     nY,minY,maxY);
     lName.str("");
     lName << "p_dRmin_" << iL;
     p_dRmin[iL] = new TH1F(lName.str().c_str(),";dRmin;events",
			    100,0,1);
   }

   p_hitMeanPuContrib = new TH2F("p_hitMeanPuContrib",";layer;E_{PU} (MIPs) from mean;hits",nLayers_,0,nLayers_,1000,0,50);
   p_hitMeanPuContrib->StatOverflows();

   p_hitEventPuContrib = new TH2F("p_hitEventPuContrib",";layer;E_{PU} (MIPs) from RC;hits",nLayers_,0,nLayers_,1000,0,50);
   p_hitEventPuContrib->StatOverflows();

   p_diffPuContrib = new TH1F("p_diffPuContrib",";E_{PU}^{RC}-E_{PU}^{avg} (MIPs) in SR2;events",1000,-2000,2000);
   p_diffPuContrib->StatOverflows();

   p_etavsphi = new TH2F("p_etavsphi",";#phi_{hit};#eta_{hit};n_{hits}",
			 900,-3.1416,3.1416,
			 250,1.4,3.0);

   p_etavsphi_max = new TH2F("p_etavsphi_max",";#phi_{max};#eta_{max};n_{events}",
			 900,-3.1416,3.1416,
			 250,1.4,3.0);
   p_etavsphi_truth = new TH2F("p_etavsphi_truth",";#phi_{gen};#eta_{gen};n_{events}",
			 900,-3.1416,3.1416,
			 250,1.4,3.0);

   p_yvsx_max = new TH2F("p_yvsx_max",";x_{max};y_{max};n_{events}",
			 nX,minX,maxX,nY,minY,maxY);
   p_yvsx_truth = new TH2F("p_yvsx_truth",";x_{gen};y_{gen};n_{events}",
			   nX,minX,maxX,nY,minY,maxY);


   p_residuals_x = new TH1F("p_residuals_x",";xreco-xtruth (mm)",1000,-50,50);
   p_residuals_y = new TH1F("p_residuals_y",";yreco-ytruth (mm)",1000,-50,50);
   p_residuals_x->StatOverflows();
   p_residuals_y->StatOverflows();

 }

 void PositionFit::initialiseFitHistograms(){

   outputFile_->cd(outputDir_.c_str());

   //check if already defined
   if (!p_chi2[0]){

     p_recoXvsLayer = new TH2F("p_recoXvsLayer",";layer;weighted x (mm);n_{events}",nLayers_,0,nLayers_,1500,-1500,1500);
     p_recoYvsLayer = new TH2F("p_recoYvsLayer",";layer;weighted y (mm);n_{events}",nLayers_,0,nLayers_,1500,-1500,1500);
     p_recoZvsLayer = new TH2F("p_recoZvsLayer",";layer;avg z (mm);n_{events}",nLayers_,0,nLayers_,3000,3170,3470);
     p_truthXvsLayer = new TH2F("p_truthXvsLayer",";layer;weighted x (mm);n_{events}",nLayers_,0,nLayers_,1500,-1500,1500);
     p_truthYvsLayer = new TH2F("p_truthYvsLayer",";layer;weighted x (mm);n_{events}",nLayers_,0,nLayers_,1500,-1500,1500);
     p_fitXvsLayer = new TH2F("p_fitXvsLayer",";layer;fit x (mm);n_{events}",nLayers_,0,nLayers_,1500,-1500,1500);
     p_fitYvsLayer = new TH2F("p_fitYvsLayer",";layer;fit y (mm);n_{events}",nLayers_,0,nLayers_,1500,-1500,1500);
     p_nLayersFit = new TH1F("p_nLayersFit",";#layers in fit;n_{events}",31,-0.5,30.5);
     p_iterations = new TH1F("p_iterations",";#iterations;n_{events}",100,0,100);

     p_chi2[0] = new TH1F("p_chi2",";#chi^{2};n_{events}",1000,0,2000);
     p_chi2[1] = new TH1F("p_chi2_truth",";#chi^{2};n_{events}",100,0,500);
     p_chi2overNDF[0] = new TH1F("p_chi2overNDF",";#chi^{2}/NDF;n_{events}",1000,0,100);
     p_chi2overNDF[1] = new TH1F("p_chi2overNDF_truth",";#chi^{2}/NDF;n_{events}",10,0,10);
     for (unsigned rt(0); rt<2;++rt){
       //p_chi2[rt]->StatOverflows();
       //p_chi2overNDF[rt]->StatOverflows();
     }

     p_impactX0[0] = new TH1F("p_impactX0",";x at z=0 (mm);n_{events}",500,-250,250);
     p_impactX0[1] = new TH1F("p_impactX0_truth",";x at z=0 (mm);n_{events}",500,-250,250);
     p_impactY0[0] = new TH1F("p_impactY0",";y at z=0 (mm);n_{events}",500,-250,250);
     p_impactY0[1] = new TH1F("p_impactY0_truth",";y at z=0 (mm);n_{events}",500,-250,250);
     p_impactXFF[0] = new TH1F("p_impactXFF",";x front face impact (mm);n_{events}",1500,-1500,1500);
     p_impactXFF[1] = new TH1F("p_impactXFF_truth",";x front face impact (mm);n_{events}",1500,-1500,1500);
     p_impactYFF[0] = new TH1F("p_impactYFF",";y front face impact (mm);n_{events}",1500,-1500,1500);
     p_impactYFF[1] = new TH1F("p_impactYFF_truth",";y front face impact (mm);n_{events}",1500,-1500,1500);
     p_impactX14[0] = new TH1F("p_impactX14",";x layer 14 impact (mm);n_{events}",1500,-1500,1500);
     p_impactX14[1] = new TH1F("p_impactX14_truth",";x layer 14 impact (mm);n_{events}",1500,-1500,1500);
     p_impactY14[0] = new TH1F("p_impactY14",";y layer 14 impact (mm);n_{events}",1500,-1500,1500);
     p_impactY14[1] = new TH1F("p_impactY14_truth",";y layer 14 impact (mm);n_{events}",1500,-1500,1500);
     p_impactZx[0] = new TH1F("p_impactZx",";z_{x} impact (mm);n_{events}",1000,-500,500);
     p_impactZx[1] = new TH1F("p_impactZx_truth",";z_{x} impact (mm);n_{events}",1000,-500,500);
     p_impactZy[0] = new TH1F("p_impactZy",";z_{y} impact (mm);n_{events}",1000,-500,500);
     p_impactZy[1] = new TH1F("p_impactZy_truth",";z_{y} impact (mm);n_{events}",1000,-500,500);
     p_tanAngleX[0] = new TH1F("p_tanAngleX",";x direction tanAngle (rad);n_{events}",500,-1,1);
     p_tanAngleX[1] = new TH1F("p_tanAngleX_truth",";x direction tanAngle (rad);n_{events}",500,-1,1);
     p_tanAngleY[0] = new TH1F("p_tanAngleY",";y direction tanAngle (rad);n_{events}",500,-1,1);
     p_tanAngleY[1] = new TH1F("p_tanAngleY_truth",";y direction tanAngle (rad);n_{events}",500,-1,1);

     p_impactX0_residual = new TH1F("p_impactX0_residual",";residual x z=0 impact (mm);n_{events}",200,-200,200);
     p_impactXFF_residual = new TH1F("p_impactXFF_residual",";residual x front face impact (mm);n_{events}",200,-10,10);
     p_impactX14_residual = new TH1F("p_impactX14_residual",";residual x layer 14 impact (mm);n_{events}",200,-10,10);
     p_tanAngleX_residual = new TH1F("p_tanAngleX_residual",";residual x direction tanAngle (rad);n_{events}",200,-0.1,0.1);
     p_angleX_residual = new TH1F("p_angleX_residual",";residual x direction angle (rad);n_{events}",200,-0.1,0.1);
     p_impactY0_residual = new TH1F("p_impactY0_residual",";residual y z=0 impact (mm);n_{events}",200,-200,200);
     p_impactYFF_residual = new TH1F("p_impactYFF_residual",";residual y front face impact (mm);n_{events}",200,-10,10);
     p_impactY14_residual = new TH1F("p_impactY14_residual",";residual y layer 14 impact (mm);n_{events}",200,-10,10);
     p_impactZx_residual = new TH1F("p_impactZx_residual",";residual z at (x=0) (mm);n_{events}",500,-500,500);
     p_impactZy_residual = new TH1F("p_impactZy_residual",";residual z at (y=0) (mm);n_{events}",500,-500,500);
     p_tanAngleY_residual = new TH1F("p_tanAngleY_residual",";residual y direction tanAngle (rad);n_{events}",200,-0.1,0.1);
     p_angleY_residual = new TH1F("p_angleY_residual",";residual y direction angle (rad);n_{events}",200,-0.1,0.1);

     p_eta_reco = new TH1F("p_eta_reco",";reco #eta;n_{events}",200,1.4,3.0);
     p_phi_reco = new TH1F("p_phi_reco",";reco #phi (rad);n_{events}",200,-3.1416,3.1416);
     p_eta_truth = new TH1F("p_eta_truth",";truth #eta;n_{events}",200,1.4,3.0);
     p_phi_truth = new TH1F("p_phi_truth",";truth #phi (rad);n_{events}",200,-3.1416,3.1416);
     p_eta_residual = new TH1F("p_eta_residual",";residual #eta;n_{events}",200,-0.1,0.1);
     p_phi_residual = new TH1F("p_phi_residual",";residual #phi (rad);n_{events}",200,-0.1,0.1);


     p_positionReso = new TH1F("p_positionReso",";#sigma_{x,y} (mm);n_{events}",500,0,50);
     p_angularReso = new TH1F("p_angularReso",";#sigma_{#theta} (rad);n_{events}",500,0,1);
   }
 }

 bool PositionFit::getGlobalMaximum(const unsigned ievt, 
				    const unsigned nVtx, 
				    std::vector<HGCSSRecoHit> *rechitvec,
				    const ROOT::Math::XYZVector & truthPos0, 
				    double & aPhimax,double & aEtamax){

   bool oneresult = true;
   TH2F *letavsphi = new TH2F("letavsphi",";#phi;#eta;hits",
			      900,-3.1416,3.1416,
			      250,1.4,3.0);

   /*
   TH2F *letavsphi_layer[nLayers_];
   if (nVtx>0){
     std::ostringstream lName;
     for (unsigned iL(0); iL<nLayers_; ++iL){
       lName.str("");
       lName << "letavsphi_layer" << iL;
       letavsphi_layer[iL] = new TH2F(lName.str().c_str(),";#phi;#eta;hits",
				      900,-3.1416,3.1416,
				      250,1.4,3.0);
     }
   }
   */

   for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
     const HGCSSRecoHit & lHit = (*rechitvec)[iH];

     double posx = lHit.get_x();
     if (fixForPuMixBug_) posx-=1.25;
     double posy = lHit.get_y();
     if (fixForPuMixBug_) posy-=1.25;
     double posz = lHit.get_z();
     //double radius = sqrt(posx*posx+posy*posy);
     double energy = lHit.energy();
     //unsigned layer = lHit.layer();
     //if (energy>1) std::cout << "Hit " << layer << " " << posx << " " << posy << " " << posz << " " << energy << std::endl;

     if (debug_>1) {
       std::cout << " --  RecHit " << iH << "/" << (*rechitvec).size() << " --" << std::endl
		 << " --  position x,y " << posx << "," << posy << std::endl;
       lHit.Print(std::cout);
     }

     ROOT::Math::XYZVector pos(posx,posy,posz);
     letavsphi->Fill(pos.phi(),pos.eta(),energy);

     //if (nVtx>0){
     //letavsphi_layer[layer]->Fill(pos.phi(),pos.eta(),energy);
     //}
     p_etavsphi->Fill(pos.phi(),pos.eta(),energy);
   }//loop on hits

   if (debug_)  std::cout << std::endl;

   //add histograms for all layers but iL
   //for (unsigned iL(0); iL<nLayers_; ++iL){
   //letavsphi_layer[iL]->Add(letavsphi,letavsphi_layer[iL],1,-1);
   //}

   //get position of maximum E tower
   int maxbin = -1;
   int binx=-1,biny=-1,binz=-1;
   //compare with truth
   double etatruth = truthPos0.eta();
   double phitruth = truthPos0.phi();
   double dRmin = 10;
   int binxmin=-1,binymin=-1;
   unsigned counter = 0;
   while (dRmin>maxdR_) {
     letavsphi->SetBinContent(maxbin,0);
     maxbin = letavsphi->GetMaximumBin();
     letavsphi->GetBinXYZ(maxbin,binx,biny,binz);
     double leta = letavsphi->GetYaxis()->GetBinCenter(biny);
     double lphi = letavsphi->GetXaxis()->GetBinCenter(binx);
     double deta = leta-etatruth;
     double dphi = DeltaPhi(lphi,phitruth);
     double dR = sqrt(pow(deta,2)+pow(dphi,2));
     if (dR<dRmin){
       dRmin=dR;
       binxmin=binx;
       binymin=biny;
     }
     if (counter>200) {
       break;
     }
     counter++;
   }

   if (debug_ || (counter>200)) std::cout << " -- Maximum found after " << counter << " trials, dRmin = " << dRmin;
   if (counter>200) std::cout << " --> away from truth pos for evt " << ievt ;
   if (debug_ || (counter>200)) std::cout << std::endl;  

   if (counter>1) oneresult = false;
   aPhimax =letavsphi->GetXaxis()->GetBinCenter(binxmin); 
   aEtamax =letavsphi->GetYaxis()->GetBinCenter(binymin); 

   //p_numberOfMaxTried->Fill(counter);
   //p_dRMaxTruth->Fill(dRmin);
   /*
   if (nVtx>0){
     double binxsize =  letavsphi->GetXaxis()->GetBinWidth(binx);
     double binysize =  letavsphi->GetYaxis()->GetBinWidth(biny);

     //allow for +/- 1 bin in each direction
     double dRmax = 2*sqrt(pow(binxsize,2)+pow(binysize,2));

     //unsigned counter = 0;
     //std::vector<double> layervec;
     //std::vector<double> etavec;
     //std::vector<double> phivec;
     std::vector<std::pair<double,double> > etaphivec;

     for (unsigned iL(0); iL<nLayers_; ++iL){
       int maxbin_layer = letavsphi_layer[iL]->GetMaximumBin();
       int binxl,binyl,binzl;
       letavsphi_layer[iL]->GetBinXYZ(maxbin_layer,binxl,binyl,binzl);
       double leta = letavsphi_layer[iL]->GetYaxis()->GetBinCenter(binyl);
       double lphi = letavsphi_layer[iL]->GetXaxis()->GetBinCenter(binxl);
       etaphivec.push_back(std::pair<double,double>(leta,lphi));
       //double deta = leta-aEtamax;
       //double dphi = DeltaPhi(lphi,aPhimax);
       //double dR = sqrt(pow(deta,2)+pow(dphi,2));
       //if (dR>dRmax) {
       //layervec.push_back(iL);
       //etavec.push_back(leta);
       //phivec.push_back(lphi);
       //counter++;
       //}
     }

     std::pair<unsigned, std::pair<double,double> > lmaj = findMajorityValue(etaphivec);
     double leta = lmaj.second.first;
     double lphi = lmaj.second.second;
     double deta = leta-aEtamax;
     double dphi = DeltaPhi(lphi,aPhimax);
     double dR = sqrt(pow(deta,2)+pow(dphi,2));

     if (dR>dRmax) {
       std::cout << " -- Warning ! Event " << ievt << " with " << lmaj.first << " majority layers dRmax=" << dRmax << " away, probably from PU." << std::endl
		 << " All layers 0-29 eta-phi = " << aEtamax << " " << aPhimax << std::endl
		 << " Maj value = " << leta << " " << lphi << std::endl;
       oneresult = false;
     }
     // if (counter > 0){
     //std::cout << " -- Warning ! Event " << ievt << " with " << counter << " layers dRmax=" << dRmax << " away, probably from PU." << std::endl
     // 		<< " All layers 0-29 eta-phi = " << aEtamax << " " << aPhimax << std::endl;
     //   for (unsigned iC(0); iC<counter;++iC){
     // 	std::cout << " Removing layer " << layervec[iC] << " eta-phi = " << etavec[iC] << " " << phivec[iC] << std::endl;
     //   }
   //}
   }//if nvtx>0
   */

   if (debug_) std::cout << " MaxE cell eta,phi = " << aEtamax << " " << aPhimax << std::endl;

   letavsphi->Delete();
   /*if (nVtx>0) {
     for (unsigned iL(0); iL<nLayers_; ++iL){
       letavsphi_layer[iL]->Delete();
     }
     }*/

   return oneresult;
 }

 void PositionFit::findSeeds(std::vector<HGCSSRecoHit> *rechitvec,
			     std::vector<bool> & seedable){

   for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
     const HGCSSRecoHit & lHit = (*rechitvec)[iH];
     if (lHit.energy()>seedMipThreshold_) seedable[iH] = true;
   }

 }

 unsigned PositionFit::getClusters(std::vector<HGCSSRecoHit> *rechitvec,
				   HGCSSClusterVec & output){

   const unsigned nHits = (*rechitvec).size();
   Clusterizer lClusterizer(debug_);
   std::vector<bool> rechitMask;
   rechitMask.resize(nHits,true);
   std::vector<bool> seedable;
   seedable.resize(nHits,false);

   findSeeds(rechitvec,seedable);

   lClusterizer.buildClusters(rechitvec,rechitMask,seedable,output);

   return output.size();

 }

 bool PositionFit::setTruthInfo(std::vector<HGCSSGenParticle> *genvec, const int G4TrackID){

   bool found = false;
   for (unsigned iP(0); iP<(*genvec).size(); ++iP){//loop on gen particles    
     //if ((*genvec).size()!= 1) (*genvec)[iP].Print(std::cout);
     if ((*genvec)[iP].trackID()==G4TrackID){
       found = true;
       //set direction
       truthDir_ = Direction((*genvec)[iP].px()/(*genvec)[iP].pz(),(*genvec)[iP].py()/(*genvec)[iP].pz());
       p_etavsphi_truth->Fill(truthDir_.phi(),truthDir_.eta());
 
       double zvtx = (*genvec)[iP].z()-((*genvec)[iP].y()/truthDir_.tanangle_y);
       double xvtx = (*genvec)[iP].x()-(((*genvec)[iP].z()-zvtx)*truthDir_.tanangle_x);
       truthVtx_ = ROOT::Math::XYZPoint(xvtx,0,zvtx);
       //in GeV
       truthE_ = (*genvec)[iP].E()/1000.;

       if (debug_) std::cout << " Found truth vertex pos at: x=" << xvtx << " z=" << zvtx << " xpos was: " << (*genvec)[iP].x() << std::endl;

       break;
     }

   }
   for (unsigned iL(0); iL<nLayers_; ++iL){
     p_genxy[iL]->Fill(truthPos(iL).X(),truthPos(iL).Y(),1);
   }

   if (!found){
     std::cout << " - Info: no photon G4trackID=" << G4TrackID << " found, already converted or outside of acceptance ..." << std::endl;
     //std::cout << " -- Number of genparticles: " << (*genvec).size() << std::endl;
     //for (unsigned iP(0); iP<(*genvec).size(); ++iP){//loop on gen particles 
     //std::cout << " --- particle " << iP << std::endl;
     //(*genvec)[iP].Print(std::cout);
     //}
   }
   else {
     p_nGenParticles->Fill((*genvec).size());
   }
   return found;

 }

/*
 bool PositionFit::getTruthPosition(std::vector<HGCSSGenParticle> *genvec,std::vector<ROOT::Math::XYPoint> & truthPos, const int G4TrackID){

   bool found = false;
   for (unsigned iP(0); iP<(*genvec).size(); ++iP){//loop on gen particles    
     //if ((*genvec).size()!= 1) (*genvec)[iP].Print(std::cout);
     if ((*genvec)[iP].trackID()==G4TrackID){
       if ((*genvec)[iP].pdgid() != 22) {
	 std::cout << " -- Error ! Particle trackID " << G4TrackID << " is not a photon ! pdgid = " << (*genvec)[iP].pdgid() << std::endl;
	 break;
       }
       found = true;
       double x0 = (*genvec)[iP].x();
       double y0 = (*genvec)[iP].y();
       double z0 = (*genvec)[iP].z();
       //double p = sqrt(pow((*genvec)[iP].px(),2)+pow((*genvec)[iP].py(),2)+pow((*genvec)[iP].pz(),2));
       //std::cout << "init : " << x0 << " " << y0 << " " << z0 << std::endl;
       //fill layers by propagating with momentum
       //ROOT::Math::XYZVector unit((*genvec)[iP].px()/p,(*genvec)[iP].py()/p,(*genvec)[iP].pz()/p);

       //std::cout << " Gen particle eta,phi = " << unit.eta() << " " << unit.phi() << std::endl;
       p_etavsphi_truth->Fill((*genvec)[iP].phi(),(*genvec)[iP].eta());
       //std::cout << " Truth pos eta-phi = " << (*genvec)[iP].eta() << " " << (*genvec)[iP].phi() << std::endl;

       for (unsigned iL(0); iL<nLayers_; ++iL){
	 //double xy = (avgZ_[iL]-z0)/sinh(unit.eta());
	 //double x = xy*cos(unit.phi())+x0;
	 //double y = xy*sin(unit.phi())+y0;
	 double x = x0+(avgZ_[iL]-z0)*(*genvec)[iP].px()/(*genvec)[iP].pz();
	 double y = y0+(avgZ_[iL]-z0)*(*genvec)[iP].py()/(*genvec)[iP].pz();
	 //std::cout << "Lay " << iL << ": " << x << " " << y << " " << avgZ_[iL] << std::endl;
	 p_genxy[iL]->Fill(x,y,1);
	 truthPos[iL] = ROOT::Math::XYPoint(x,y);
       }

     }
   }//loop on gen particles

   if (!found){
     std::cout << " - Info: no photon G4trackID=" << G4TrackID << " found, already converted or outside of acceptance ..." << std::endl;
     //std::cout << " -- Number of genparticles: " << (*genvec).size() << std::endl;
     //for (unsigned iP(0); iP<(*genvec).size(); ++iP){//loop on gen particles 
     //std::cout << " --- particle " << iP << std::endl;
     //(*genvec)[iP].Print(std::cout);
     //}
   }
   else {
     p_nGenParticles->Fill((*genvec).size());
   }
   return found;

 }
*/

 bool PositionFit::getZpositions(const unsigned versionNumber){
   std::ifstream fin;
   std::ostringstream finname;
   //finname << outFolder_ << "/zPositions.dat";
   finname << "data/zPositions_v" << versionNumber << ".dat";
   fin.open(finname.str());
   if (!fin.is_open()){
     std::cout << " Cannot open input file " << finname.str() << "! Refilling now..." << std::endl;
     return false;
   }

   std::cout << " Reading z position per layer from input file " << finname.str() << std::endl;

   std::vector<unsigned> layerId;
   std::vector<double> posz;
   layerId.reserve(nLayers_);
   posz.reserve(nLayers_);

   while (!fin.eof()){
     unsigned l=nLayers_;
     double z=0;
     fin>>l>>z;
     if (l<nLayers_){
       avgZ_.push_back(z);
       std::cout << " Layer " << l << ", z = " << z << std::endl;
     }
   }

   if (avgZ_.size() != nLayers_) {
     std::cout << " -- Warning! Problem in extracting z positions, did not find one value per layer. Please check input file: " << finname.str() << std::endl
	       << " Proceeding to refilling them from simtree." << std::endl;
     return false;
   }

   fin.close();
   return true;

 }

void PositionFit::getZpositions(const unsigned versionNumber,
				TTree *aSimTree,
				 const unsigned nEvts){


   std::vector<HGCSSSimHit> * simhitvec = 0;
   aSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);

   std::ofstream fout;
   std::ostringstream foutname;
   //foutname << outFolder_ << "/zPositions.dat";
   foutname << "data/zPositions_v" << versionNumber << ".dat";
   fout.open(foutname.str());
   if (!fout.is_open()){
     std::cout << " Cannot open outfile " << foutname.str() << " for writing ! Exiting..." << std::endl;
     exit(1);
   }

   std::cout << "--- Filling z positions:" << std::endl
	     << "- Processing = " << nEvts  << " events out of " << aSimTree->GetEntries() << std::endl;

   avgZ_.resize(nLayers_,0);


   for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
     if (debug_) std::cout << "... Processing entry: " << ievt << std::endl;
     else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
     aSimTree->GetEntry(ievt);
     for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on rechits
       const HGCSSSimHit & lHit = (*simhitvec)[iH];
       unsigned layer = lHit.layer();
       if (layer >= nLayers_) {
	 continue;
       }

       //discard some si layers...
       if (lHit.silayer() >= nSiLayers_) continue; 

       double posz = lHit.get_z();
       //get z position of hits
       if (avgZ_[layer]<posz) avgZ_[layer]=posz;
     }
   }

   std::cout << " --- Z positions of layers: " << std::endl;
   for (unsigned iL(0); iL<nLayers_;++iL){
     std::cout << " Layer " << iL << ", z = " << avgZ_[iL] << std::endl;
     fout << iL << " " << avgZ_[iL] << std::endl;
   }

   fout.close();

 }

 void PositionFit::getInitialPositions(TTree *aSimTree, 
				       TTree *aRecTree,
				       const unsigned nEvts,
				       const unsigned G4TrackID){

   initialiseClusterHistograms();
   initialisePositionHistograms();
   //HGCSSDetector & myDetector = theDetector();

   //////////////////////////////////////////////////
   ///////// Event loop /////////////////////////////
   //////////////////////////////////////////////////

   HGCSSEvent * event = 0;
   std::vector<HGCSSSamplingSection> * ssvec = 0;
   std::vector<HGCSSSimHit> * simhitvec = 0;
   std::vector<HGCSSRecoHit> * rechitvec = 0;
   std::vector<HGCSSGenParticle> * genvec = 0;
   unsigned nPuVtx = 0;

   aSimTree->SetBranchAddress("HGCSSEvent",&event);
   aSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
   aSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
   aSimTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);

   aRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
   if (aRecTree->GetBranch("nPuVtx")) aRecTree->SetBranchAddress("nPuVtx",&nPuVtx);

   std::cout << "- Processing = " << nEvts  << " events out of " << aSimTree->GetEntries() << std::endl;

   //bool firstEvent = true;

   unsigned nConvertedPhotons = 0;
   unsigned nTooFar = 0;

   //output tree
   TTree *outtree = 0;
   std::vector<double> truthPosX;
   std::vector<double> truthPosY;
   std::vector<std::vector<double> > Exy;
   std::vector<double> init;
   init.resize(25,0);
   Exy.resize(nLayers_,init);
   if (saveEtree_){
     outputFile_->cd();
     outtree = new TTree("EcellsSR2","Tree to save energies in each cell of SR2");
     truthPosX.resize(nLayers_,0);
     truthPosY.resize(nLayers_,0);
     for (unsigned iL(0);iL<nLayers_;++iL){
       std::ostringstream label;
       label.str("");     
       label << "TruthPosX_" << iL;
       outtree->Branch(label.str().c_str(),&truthPosX[iL]);
       label.str("");     
       label << "TruthPosY_" << iL;
       outtree->Branch(label.str().c_str(),&truthPosY[iL]);
       for (unsigned iy(0);iy<5;++iy){
	 for (unsigned ix(0);ix<5;++ix){
	   unsigned idx = 5*iy+ix;
	   label.str("");     
	   label << "E_" << iL << "_" << idx;
	   outtree->Branch(label.str().c_str(),&Exy[iL][idx]);
	 }
       }
     }
   }


  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug_) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
    
    aSimTree->GetEntry(ievt);

    aRecTree->GetEntry(ievt);

    if (debug_) std::cout << " nPuVtx = " << nPuVtx << std::endl;

    if (debug_){
      std::cout << "... Size of hit vectors: sim = " <<  (*simhitvec).size() << ", reco = " << (*rechitvec).size()<< " gen " << (*genvec).size() << std::endl;
    }

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    //////// output files to save position for chi2 fit //////////
    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////

    bool found = setTruthInfo(genvec,G4TrackID);

    //std::vector<ROOT::Math::XYPoint> truthPos;
    //truthPos.resize(nLayers_,ROOT::Math::XYPoint(0,0));
    //bool found = getTruthPosition(genvec,truthPos,G4TrackID);
    if (!found) {
      nConvertedPhotons++;
      continue;
    }
    if (saveEtree_) {
      for (unsigned iL(0);iL<nLayers_;++iL){
	truthPosX[iL] = truthPos(iL).X();
	truthPosY[iL] = truthPos(iL).Y();
      }
    }
    
    p_yvsx_truth->Fill(truthPos(10).X(),truthPos(10).Y());


    //get position of maximum to start with

    //from global max -- not working in PU
    /*
    double phimax = 0;
    double etamax = 0;
    ROOT::Math::XYZVector truthPos0(truthPos(0).X(),truthPos(0).Y(),avgZ_[0]);
    bool oneresult = getGlobalMaximum(ievt,nPuVtx,rechitvec,truthPos0,phimax,etamax);
    if (!oneresult) nMultipleMax++;
    */

    if (!getInitialPosition(ievt,nPuVtx,rechitvec,nTooFar,Exy)) continue;
    
    if (saveEtree_) outtree->Fill();

    //firstEvent = false;
  }//loop on entries

  std::cout << " -- Number of converted photons: " << nConvertedPhotons << std::endl;
  std::cout << " -- Number of events with closest cluster away from truth within dR " << maxdR_ << " : " << nTooFar << std::endl;
  
 }

bool PositionFit::getInitialPosition(const unsigned ievt,
				     const unsigned nPuVtx, 
				     std::vector<HGCSSRecoHit> *rechitvec,
				     unsigned & nTooFar,
				     std::vector<std::vector<double> > & Exy){
  
  //from clusters -- using Lindsey's clustering
  HGCSSClusterVec lClusVec;
  unsigned nClusters = getClusters(rechitvec,lClusVec);
  
  p_nClusters->Fill(nClusters);
  
  if (nClusters == 0) {
    std::cout << " -- Error, no cluster found !" << std::endl;
    exit(1);
    //continue;
  }
  else if (debug_) std::cout << " -- Found " << nClusters << " clusters." << std::endl;
  
  //loop over clusters to find closest to the truth
  //get eta-phi of the seed ...
  double etatruth = truthDir_.eta();
  double phitruth = truthDir_.phi();
  double dRmin = 10;
  unsigned clusIdx = 0;
  for (unsigned iClus(0); iClus<nClusters;++iClus){
    const HGCSSCluster & lCluster = lClusVec[iClus];
    double leta = lCluster.eta();
    double lphi = lCluster.phi();
    p_clusnHits_all->Fill(lCluster.nRecHits());
    if (lCluster.energy()>0) p_seedEoverE_all->Fill(lCluster.getSeedE()/lCluster.energy());
    p_clusLayer_all->Fill(lCluster.layer());
    p_clusWidth_all->Fill(lCluster.width());
    p_etavsphi->Fill(lphi,leta,lCluster.energy());
    double deta = leta-etatruth;
    double dphi = DeltaPhi(lphi,phitruth);
    double dR = sqrt(pow(deta,2)+pow(dphi,2));
    if (dR<dRmin){
      dRmin = dR;
      clusIdx = iClus;
    }
    double detas = leta-lCluster.getSeedEta();
    double dphis = DeltaPhi(lphi,lCluster.getSeedPhi());
    p_seeddeta_all->Fill(detas);
    p_seeddphi_all->Fill(dphis);
    if (debug_) {
      std::cout << " Cluster " << iClus << ":" ;
      lCluster.Print(std::cout);
    }
  }
  
  const HGCSSCluster & lCluster = lClusVec[clusIdx];
  double phimax = lCluster.phi();//getSeedPhi();
  double etamax = lCluster.eta();//getSeedEta();
  
  if (debug_){
    std::cout << " -- evt " << ievt << " found cluster mindR= " << dRmin;
    lCluster.Print(std::cout);
    std::cout << " Truth pos eta-phi = " << etatruth << " " << phitruth << std::endl;
  }
  
  p_mindRtruth->Fill(dRmin);
  
  if (dRmin > 0.1) {
    nTooFar++;
    return false;
  }
  
  //fill selected cluster histos
  p_clusnHits_sel->Fill(lCluster.nRecHits());
  if (lCluster.energy()>0) p_seedEoverE_sel->Fill(lCluster.getSeedE()/lCluster.energy());
  p_clusLayer_sel->Fill(lCluster.layer());
  p_clusWidth_sel->Fill(lCluster.width());
  double detas = lCluster.eta()-lCluster.getSeedEta();
  double dphis = DeltaPhi(lCluster.phi(),lCluster.getSeedPhi());
  p_seeddeta_sel->Fill(detas);
  p_seeddphi_sel->Fill(dphis);
  p_etavsphi_max->Fill(phimax,etamax);
  
  std::vector<double> xmax;
  xmax.resize(nLayers_,0);
  std::vector<double> ymax;
  ymax.resize(nLayers_,0);
  //getMaximumCellFromGeom(phimax,etamax,xmax,ymax);
  getMaximumCell(rechitvec,phimax,etamax,xmax,ymax);
  p_yvsx_max->Fill(xmax[10],ymax[10]);
  
  //get PU contrib from elsewhere in the event
  //loop over phi with same etamax
  //take average per layer: not all 9 cells of 3*3 area have hits...
  std::vector<double> puE;
  puE.resize(nLayers_,0);
  if (nPuVtx>0){
    unsigned nRandomCones = 50;
    double phistep = TMath::Pi()/nRandomCones;
    if (debug_) std::cout << "--- etamax = " << etamax << " phimax=" << phimax << " phistep = " << phistep << std::endl;
    for (unsigned ipm(0);ipm<nRandomCones;++ipm){
      std::vector<double> xmaxrc;
      xmaxrc.resize(nLayers_,0);
      std::vector<double> ymaxrc;
      ymaxrc.resize(nLayers_,0);
      double phirc = phimax-TMath::Pi();
      if (phirc < -1.*TMath::Pi()) phirc+=2.*TMath::Pi();
      if (phirc > TMath::Pi()) phirc-=2.*TMath::Pi();
      if (ipm%2==0) phirc += ipm/2*phistep+phistep/2.;
      else  phirc = phirc - ipm/2*phistep-phistep/2.;
      if (phirc < -1.*TMath::Pi()) phirc+=2.*TMath::Pi();
      if (phirc > TMath::Pi()) phirc-=2.*TMath::Pi();
      //take from geom to not be biased by hit having PU, because
      //not from geom means find cell with a hit closest to maxpos...
      getMaximumCellFromGeom(phirc,etamax,xmaxrc,ymaxrc);
      if (debug_>1) std::cout << "rc #" << ipm << " phirc=" << phirc << " xmax[10]=" << xmaxrc[10] << " ymax[10]=" << ymaxrc[10] << " r=" << sqrt(pow(xmaxrc[10],2)+pow(ymaxrc[10],2)) << std::endl;
      getPuContribution(rechitvec,xmaxrc,ymaxrc,puE);
    }
    
    //normalise to one cell: must count cells with 0 hit !
    //use cell size at etamax...
    for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
      if (debug_) std::cout << "layer " << iL ;
      unsigned nCells = nRandomCones*pow(nSR_*geomConv_.cellSize()/geomConv_.cellSize(iL,etamax),2);
      puE[iL] = puE[iL]/nCells;
      if (debug_) std::cout << " Epu=" << puE[iL] << std::endl;	
    }
    
  }//if PU
  
  
  std::vector<ROOT::Math::XYPoint> recoPos;
  std::vector<double> recoE;
  recoPos.resize(nLayers_,ROOT::Math::XYPoint(0,0));
  recoE.resize(nLayers_,0);
  std::vector<unsigned> nHits;
  nHits.resize(nLayers_,0);
  
  //get energy-weighted position and energy around maximum
  getEnergyWeightedPosition(rechitvec,nPuVtx,xmax,ymax,recoPos,recoE,nHits,puE,Exy);
  
  std::ofstream fout;
  std::ostringstream foutname;
  foutname << outFolder_ << "/initialPos_evt" << ievt << ".dat";
  fout.open(foutname.str());
  if (!fout.is_open()){
    std::cout << " Cannot open outfile " << foutname.str() << " for writing ! Exiting..." << std::endl;
    exit(1);
  }
  
  if (debug_) std::cout << " Summary of reco and truth positions:" << std::endl;
  for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
    if (debug_) std::cout << iL << " nHits=" << nHits[iL] << " Max=(" << xmax[iL] << "," << ymax[iL] << ")\t Reco=(" << recoPos[iL].X() << "," << recoPos[iL].Y() << ")\t Truth=(" << truthPos(iL).X() << "," << truthPos(iL).Y() << ")" << std::endl;
    fout << iL << " " << recoPos[iL].X() << " " << recoPos[iL].Y() << " " << truthPos(iL).X() << " " << truthPos(iL).Y() ;
    if (!doMatrix_) fout << " " << recoE[iL];
    fout << std::endl;
  }

  fout.close();

  if (doMatrix_) fillErrorMatrix(recoPos,nHits);

  return true;

}

void PositionFit::getMaximumCellFromGeom(const double & phimax,const double & etamax,std::vector<double> & xmax,std::vector<double> & ymax){

  for (unsigned iL(0); iL<nLayers_;++iL){
    double theta = 2*atan(exp(-1.*etamax));
    double rho = avgZ_[iL]/cos(theta);
    xmax[iL] = rho*sin(theta)*cos(phimax);
    ymax[iL] = rho*sin(theta)*sin(phimax);

    if (xmax[iL]>0) xmax[iL]=static_cast<int>((xmax[iL]+4.999999)/10.)*10;
    else xmax[iL]=static_cast<int>((xmax[iL]-4.999999)/10.)*10;
    if (ymax[iL]>0) ymax[iL]=static_cast<int>((ymax[iL]+4.999999)/10.)*10;
    else ymax[iL]=static_cast<int>((ymax[iL]-4.999999)/10.)*10;

  }//loop on layers

}

void PositionFit::getMaximumCell(std::vector<HGCSSRecoHit> *rechitvec,const double & phimax,const double & etamax,std::vector<double> & xmax,std::vector<double> & ymax){
  
  //std::vector<double> dRmin;
  //dRmin.resize(nLayers_,10);
  std::vector<double> xmaxgeom;
  xmaxgeom.resize(nLayers_,0);
  std::vector<double> ymaxgeom;
  ymaxgeom.resize(nLayers_,0);
  getMaximumCellFromGeom(phimax,etamax,xmaxgeom,ymaxgeom);

  //choose cell with maximum energy from 3*3 array around geom pos of max.
  std::vector<double> Emax;
  Emax.resize(nLayers_,0);

  for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on rechits
    const HGCSSRecoHit & lHit = (*rechitvec)[iH];
    
    unsigned layer = lHit.layer();
    if (layer >= nLayers_) {
      continue;
    }
    
    double posx = lHit.get_x();
    if (fixForPuMixBug_) posx-=1.25;
    double posy = lHit.get_y();
    if (fixForPuMixBug_) posy-=1.25;
    //double posz = avgZ_[layer];//lHit.get_z();
    if (debug_>1) {
      std::cout << " --  RecoHit " << iH << "/" << (*rechitvec).size() << " --" << std::endl
		<< " --  position x,y " << posx << "," << posy << std::endl;
      lHit.Print(std::cout);
    }
    double energy = lHit.energy();

    if (fabs(posx-xmaxgeom[layer]) <= 15 && 
	fabs(posy-ymaxgeom[layer]) <= 15){

      if (energy>Emax[layer]){
	Emax[layer] = energy;
	xmax[layer] = posx;
	ymax[layer] = posy;
      }
    }
      //ROOT::Math::XYZVector pos(posx,posy,posz);
      //double deta = fabs(pos.eta()-etamax);
      //double dphi = DeltaPhi(pos.phi(),phimax);
    
      //double dR = sqrt(pow(deta,2)+pow(dphi,2));
      //if (dR<dRmin[layer]) {
      //dRmin[layer] = dR;
      //xmax[layer] = posx;
      //ymax[layer] = posy;
      //}
    
    p_recoxy[layer]->Fill(posx,posy,energy);
    
  }//loop on rechits
    
  //for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
  //p_dRmin[iL]->Fill(dRmin[iL]);
  //}
}

void PositionFit::getEnergyWeightedPosition(std::vector<HGCSSRecoHit> *rechitvec,
					    const unsigned nPU, 
					    const std::vector<double> & xmax,
					    const std::vector<double> & ymax,
					    std::vector<ROOT::Math::XYPoint> & recoPos,
					    std::vector<double> & recoE,
					    std::vector<unsigned> & nHits,
					    std::vector<double> & puE,
					    std::vector<std::vector<double> > & Exy,
					    const bool puSubtracted){
  
  double eSum_puMean = 0;
  double eSum_puEvt = 0;

  double steplarge = geomConv_.cellSize()*20/2.+0.1;//+0.1 to accomodate double precision
  double step = geomConv_.cellSize()*nSR_/2.+0.1;//+0.1 to accomodate double precision
  if (debug_) std::cout << "step = " << step << std::endl;

  for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on rechits
    const HGCSSRecoHit & lHit = (*rechitvec)[iH];
    double energy = lHit.energy();//in MIP already...
    unsigned layer = lHit.layer();
    double posx = lHit.get_x();
    if (fixForPuMixBug_) posx-=1.25;
    double posy = lHit.get_y();
    if (fixForPuMixBug_) posy-=1.25;


    if (fabs(posx-xmax[layer]) < steplarge && 
	fabs(posy-ymax[layer]) < steplarge){
      if (puSubtracted) {
	double leta = lHit.eta();
	if (debug_>1) std::cout << " -- Hit " << iH << ", eta=" << leta << ", energy before PU subtraction: " << energy << " after: " ;
	double lCorMean =  puDensity_.getDensity(leta,layer,geomConv_.cellSizeInCm(layer,leta),nPU);
	eSum_puMean += lCorMean;
	p_hitMeanPuContrib->Fill(layer,lCorMean);
	double lCorEvent = puE[layer];
	eSum_puEvt += lCorEvent;
	p_hitEventPuContrib->Fill(layer,lCorEvent);
	double lCor = 0;
	if (useMeanPU_) lCor = lCorMean;
	else lCor = lCorEvent;
	energy = std::max(0.,energy - lCor);
	if (debug_>1) std::cout << energy << std::endl;
      }
      //fill anyway linear weight, decide in next loop which to use
      if (fabs(posx-xmax[layer]) < step && 
	  fabs(posy-ymax[layer]) < step){
	recoPos[layer].SetX(recoPos[layer].X() + posx*energy);
	recoPos[layer].SetY(recoPos[layer].Y() + posy*energy);
	recoE[layer] += energy;
	if (energy>0) nHits[layer]++;
      }
      
      if (doLogWeight_){
	int ix = (posx-xmax[layer])/10.;
	int iy = (posy-ymax[layer])/10.;
	unsigned idx = 0;
	if ((ix > 2 || ix < -2) || (iy>2 || iy<-2)) {
	  std::cout << " error, check ix=" << ix << " iy=" << iy << " posx,y-max=" << posx-xmax[layer] << " " << posy-ymax[layer] << " step " << step << std::endl;
	  continue;
	}
	else 
	  idx = 5*(iy+2)+(ix+2);
	Exy[layer][idx] = energy;
      }
    }
    
  }//loop on rechits

  p_diffPuContrib->Fill(eSum_puEvt-eSum_puMean);

  for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
    if (nHits[iL]==0) continue;

    if (doLogWeight_){
      double Etot = 0;
      double Ex[5] = {0,0,0,0,0};
      double Ey[5] = {0,0,0,0,0};
      for (unsigned idx(0);idx<25;++idx){
	if (iL>22) Etot += Exy[iL][idx];
	else if ((idx>5 && idx<9)||
		 (idx>10 && idx<14)||
		 (idx>15 && idx<19)) Etot += Exy[iL][idx];
      }
      if (iL>22){
	Ex[0] = Exy[iL][0]+Exy[iL][5]+Exy[iL][10]+Exy[iL][15]+Exy[iL][20];
	Ex[1] = Exy[iL][1]+Exy[iL][6]+Exy[iL][11]+Exy[iL][16]+Exy[iL][21];
	Ex[2] = Exy[iL][2]+Exy[iL][7]+Exy[iL][12]+Exy[iL][17]+Exy[iL][22];
	Ex[3] = Exy[iL][3]+Exy[iL][8]+Exy[iL][13]+Exy[iL][18]+Exy[iL][23];
	Ex[4] = Exy[iL][4]+Exy[iL][9]+Exy[iL][14]+Exy[iL][19]+Exy[iL][24];
	Ey[0] = Exy[iL][0]+Exy[iL][1]+Exy[iL][2]+Exy[iL][3]+Exy[iL][4];
	Ey[1] = Exy[iL][5]+Exy[iL][6]+Exy[iL][7]+Exy[iL][8]+Exy[iL][9];
	Ey[2] = Exy[iL][10]+Exy[iL][11]+Exy[iL][12]+Exy[iL][13]+Exy[iL][14];
	Ey[3] = Exy[iL][15]+Exy[iL][16]+Exy[iL][17]+Exy[iL][18]+Exy[iL][19];
	Ey[4] = Exy[iL][20]+Exy[iL][21]+Exy[iL][22]+Exy[iL][23]+Exy[iL][24];
      }
      else {
	Ex[0] = Exy[iL][6]+Exy[iL][11]+Exy[iL][16];
	Ex[1] = Exy[iL][7]+Exy[iL][12]+Exy[iL][17];
	Ex[2] = Exy[iL][8]+Exy[iL][13]+Exy[iL][18];
	Ey[0] = Exy[iL][6]+Exy[iL][7]+Exy[iL][8];
	Ey[1] = Exy[iL][11]+Exy[iL][12]+Exy[iL][13];
	Ey[2] = Exy[iL][16]+Exy[iL][17]+Exy[iL][18];
      }
      double wx[6];
      double wy[6];
      for (unsigned i(0);i<6;++i){
	wx[i] = 0;
	wy[i] = 0;
      }
      double w0 = getW0(iL);
      //std::cout << iL << "&" << w0 << "\\\\ \n";
      for (unsigned i(0);i<5;++i){
	wx[i] = std::max(0.,log(Ex[i]/Etot)+w0);
	wy[i] = std::max(0.,log(Ey[i]/Etot)+w0);
	wx[5] += wx[i];
	wy[5] += wy[i];
      }
      double x = xmax[iL];
      //if none pass, discard layer
      if (wx[5]!=0) {
	if (iL>22) x += 10*(2*wx[4]+wx[3]-wx[1]-2*wx[0])/wx[5];
	else x += 10*(wx[2]-wx[0])/wx[5];
      }
      else nHits[iL]=0;
      double y = ymax[iL];
      if (wy[5]!=0) {
	if (iL>22) y += 10*(2*wy[4]+wy[3]-wy[1]-2*wy[0])/wy[5];
	else y += 10*(wy[2]-wy[0])/wy[5];
      }
      else nHits[iL]=0;
      
      if (nHits[iL]!=0){
	recoPos[iL].SetX(x);
	recoPos[iL].SetY(y);
      }
    }
    else {
      recoPos[iL].SetX(recoPos[iL].X()/recoE[iL]);
      recoPos[iL].SetY(recoPos[iL].Y()/recoE[iL]);
    }
  }//loop on layers
  
}

void PositionFit::getPuContribution(std::vector<HGCSSRecoHit> *rechitvec, const std::vector<double> & xmax,const std::vector<double> & ymax,std::vector<double> & puE){

  double step = geomConv_.cellSize()*nSR_/2.+0.1;

  for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on rechits
    const HGCSSRecoHit & lHit = (*rechitvec)[iH];
    double energy = lHit.energy();//in MIP already...
    unsigned layer = lHit.layer();
    double posx = lHit.get_x();
    if (fixForPuMixBug_) posx-=1.25;
    double posy = lHit.get_y();
    if (fixForPuMixBug_) posy-=1.25;

    //std::cout << "- iH " << iH << " x=" << posx << " xmax=" << xmax[layer] << " y=" << posy << " ymax=" << ymax[layer] << " step " << step << std::endl;
    if (fabs(posx-xmax[layer]) < step && 
	fabs(posy-ymax[layer]) < step){
      if (debug_>1) std::cout << "- iH " << iH 
			     << " x=" << posx 
			     << " xmax=" << xmax[layer] 
			     << " y=" << posy 
			     << " ymax=" << ymax[layer] 
			     << " step " << step 
			     << " --- Pass, layer" << layer 
			     << " E="<< energy << std::endl;
      puE[layer] += energy;
    }

  }//loop on rechits
  //exit(1);
}

void PositionFit::fillErrorMatrix(const std::vector<ROOT::Math::XYPoint> & recoPos,const std::vector<unsigned> & nHits){

  for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
    if (nHits[iL]==0) continue;
    double residual_xi = recoPos[iL].X()-truthPos(iL).X();
    double residual_yi = recoPos[iL].Y()-truthPos(iL).Y();
   if (fabs(residual_xi)>residualMax_ || fabs(residual_yi)>residualMax_) continue;
    p_residuals_x->Fill(residual_xi);
    p_residuals_y->Fill(residual_yi);
    //unsigned posmm = static_cast<unsigned>(fabs(truthPos(iL).Y())+5);
    //bool isEdge = posmm%10 <= 2 || posmm%10 >= 8;
    //if (!isEdge) continue;
    mean_[iL] += residual_xi;
    mean_[nLayers_+iL] += residual_yi;
    ++nL_mean_[iL];
    ++nL_mean_[nLayers_+iL];
    if (debug_>1) {
      std::cout << " -- Means for layer: " << iL << " " 
		<< mean_[iL] << " " 
		<< mean_[nLayers_+iL] << " "
		<< nL_mean_[iL]
		<< std::endl;
    }
     for (unsigned jL(0);jL<nLayers_;++jL){//loop on layers
      if (nHits[jL]==0) continue;
      double residual_xj = recoPos[jL].X()-truthPos(jL).X();
      double residual_yj = recoPos[jL].Y()-truthPos(jL).Y();
      if (fabs(residual_xj)>residualMax_ || fabs(residual_yj)>residualMax_) continue;
      //posmm = static_cast<unsigned>(fabs(truthPos(jL).Y())+5);
      //isEdge = posmm%10 <= 2 || posmm%10 >= 8;
      //if (!isEdge) continue;
      double sigma_x = residual_xi*residual_xj;
      double sigma_y = residual_yi*residual_yj;
      sigma_[iL][jL] += sigma_x;
      sigma_[nLayers_+iL][nLayers_+jL] += sigma_y;
      double sigma_xy = residual_xi*residual_yj;
      sigma_[iL][nLayers_+jL] += sigma_xy;
      sigma_xy = residual_xj*residual_yi;
      sigma_[nLayers_+iL][jL] += sigma_xy;

      ++nL_sigma_[iL][jL];
      ++nL_sigma_[nLayers_+iL][nLayers_+jL];
      ++nL_sigma_[iL][nLayers_+jL];
      ++nL_sigma_[nLayers_+iL][jL];
    }//loop on layers
  }//loop on layers

}


void PositionFit::finaliseErrorMatrix(){
  //finalise error matrix
  //x
  //finaliseErrorMatrix(true);
  //y
  //finaliseErrorMatrix(false);
  finaliseErrorMatrix(true);
}

void PositionFit::finaliseErrorMatrix(const bool doX){
  //finalise error matrix

  std::ofstream fmatrix;
  std::ostringstream fmatrixname;
  fmatrixname << matrixFolder_ << "/errorMatrix";
  //if (doX) fmatrixname << "_x";
  //else fmatrixname << "_y";
  fmatrixname << "_xy.dat";
  fmatrix.open(fmatrixname.str());
  if (!fmatrix.is_open()){
    std::cout << " Cannot open outfile " << fmatrixname.str() << " for writing ! Exiting..." << std::endl;
    exit(1);
  }

  //const unsigned index = (doX)? 0 : 1;

  matrix_.ResizeTo(2*nLayers_,2*nLayers_);


  //set mean values first
  for (unsigned iL(0);iL<2*nLayers_;++iL){//loop on layers
    if (nL_mean_[iL]>0) mean_[iL] = mean_[iL]/nL_mean_[iL];
    else mean_[iL] = 0;
    if (debug_>1) {
      std::cout << " -- Means for layer: " << iL << " " 
		<< mean_[iL] << " " 
		<< nL_mean_[iL]
		<< std::endl;
    }
   }
  //set sigmas and fill matrix
  for (unsigned iL(0);iL<2*nLayers_;++iL){//loop on layers
    for (unsigned jL(0);jL<2*nLayers_;++jL){//loop on layers
      if (nL_sigma_[iL][jL]>0) sigma_[iL][jL] = sigma_[iL][jL]/nL_sigma_[iL][jL];
      else sigma_[iL][jL] = 0;
     matrix_[iL][jL] = sigma_[iL][jL]-mean_[iL]*mean_[jL];
     if (debug_>1) {
       std::cout << " -- Sigma for layer: " 
		 << iL << " " << jL << " "
		 << sigma_[iL][jL] << " " 
		 << nL_sigma_[iL][jL] 
		 << " matrix = "
		 << matrix_[iL][jL]
		 << std::endl;
     }
     if (matrix_[iL][jL]!=matrix_[iL][jL]) matrix_[iL][jL] = 0;
     fmatrix << iL << " " << jL << " " << std::setprecision(15) << matrix_[iL][jL] << std::endl;

      //if (iL!=jL){
      //p_matrix->Fill(jL,iL,matrix_[index][iL][jL]);
      //}
    }
  }

  std::cout << " -- End of filling matrix" << std::endl;

  fmatrix.close();

}

void PositionFit::fillCorrelationMatrix(){
  std::cout << " -- Filling correlation matrix" << std::endl;
  outputFile_->cd(outputDir_.c_str());

  p_errorMatrix_xy = new TH2D("p_errorMatrix_xy",";i;j;M_{ij}",
			     2*nLayers_,0,2*nLayers_,
			     2*nLayers_,0,2*nLayers_);
  p_corrMatrix_xy = new TH2D("p_corrMatrix_xy",";i;j;M_{ij}",
			    2*nLayers_,0,2*nLayers_,
			    2*nLayers_,0,2*nLayers_);
  //p_errorMatrix_y = new TH2D("p_errorMatrix_y",";i;j;M_{ij}",
  //nLayers_,0,nLayers_,
  //nLayers_,0,nLayers_);
  //p_corrMatrix_y = new TH2D("p_corrMatrix_y",";i;j;M_{ij}",
  //nLayers_,0,nLayers_,
  //nLayers_,0,nLayers_);

  corrMatrix_.ResizeTo(2*nLayers_,2*nLayers_);
  for (unsigned iL(0);iL<2*nLayers_;++iL){//loop on layers
    for (unsigned jL(0);jL<2*nLayers_;++jL){//loop on layers
      p_errorMatrix_xy->Fill(iL,jL,matrix_[iL][jL]);
      if (matrix_[iL][iL]!=0 && matrix_[jL][jL]!= 0){
	corrMatrix_[iL][jL] =matrix_[iL][jL]/sqrt(matrix_[iL][iL]*matrix_[jL][jL]); 
	p_corrMatrix_xy->Fill(iL,jL,corrMatrix_[iL][jL]);
      }
    }
  }
}


bool PositionFit::fillMatrixFromFile(const bool old){
  return (fillMatrixFromFile(true,old));// && fillMatrixFromFile(false,old));
}

bool PositionFit::fillMatrixFromFile(const bool doX, const bool old){

  std::ifstream fmatrix;
  std::ostringstream fmatrixname;
  fmatrixname << matrixFolder_ << "/errorMatrix";
  if (!old){
    fmatrixname << "_xy";
    //if (doX) fmatrixname << "_x";
    //else fmatrixname << "_y";
  }
  fmatrixname << ".dat";
  fmatrix.open(fmatrixname.str());
  if (!fmatrix.is_open()){
    std::cout << " -- Cannot open outfile " << fmatrixname.str() << "! Returning false..." << std::endl;
    return false;
  }

  //const unsigned index = (doX)? 0 : 1;

  matrix_.ResizeTo(2*nLayers_,2*nLayers_);
  if (debug_) std::cout << " -- Error matrix: " << std::endl;
  while (!fmatrix.eof()){
    unsigned iL=2*nLayers_;
    unsigned jL=2*nLayers_;
    double m=0;
    fmatrix>>iL>>jL>>m;
    if (iL<2*nLayers_ && jL<2*nLayers_){
      if (debug_) std::cout << std::setprecision(15) << iL << " " << jL << " " << m << std::endl;
      matrix_[iL][jL] = m;
    }
    else if (debug_) std::cout << "!! out of bounds!" << iL << " " << jL << " " << m << std::endl;
  }
  
  if (debug_) std::cout << " -- Matrix read from file successfully." << std::endl;
  fmatrix.close();

  return true;
}

bool PositionFit::initialiseLeastSquareFit(){
  initialiseFitHistograms();

  //try reading matrix from file, if fail
  //return false and refill the matrix.
  if (doMatrix_ && !fillMatrixFromFile()) return false;
  
  //fill matrices
  if (doMatrix_) fillCorrelationMatrix();

  //get back data for each event and perform chi2 fit:
  std::cout << " -- Performing chi2 fit for each event" << std::endl;

  nInvalidFits_=0;
  nFailedFitsAfterCut_=0;
  countFailedSimFits_ = 0;
  //open new file to save accurate positions
  std::ostringstream foutname;
  foutname << outFolder_ << "/accuratePos.dat";
  fout_.open(foutname.str());
  if (!fout_.is_open()){
    std::cout << " Cannot open outfile " << foutname.str() << " for writing ! Exiting..." << std::endl;
    exit(1);
  }

  return true;

}

unsigned PositionFit::performLeastSquareFit(const unsigned ievt,
					    FitResult & fit){

    //cut outliers
  unsigned fitres = fitEvent(ievt,fit,true);
    if (fitres==1){
      // std::cout << " -- Number of genparticles: " << (*genvec).size() << std::endl;
      // for (unsigned iP(0); iP<(*genvec).size(); ++iP){//loop on gen particles 
      // 	std::cout << " --- particle " << iP << std::endl;
      // 	(*genvec)[iP].Print(std::cout);
      // }
      std::cout << " -- Event " << ievt << " skipped." << std::endl;
    }
    if (fitres>1) {
      //std::cout << " ---- Fit failed ! Reprocess event " << ievt << " cutting outliers" << std::endl;
      //nFailedFits++;
      //if (fitEvent(ievt,fit,true)>1)
      std::cout << " -- Event " << ievt << " failed fit." << std::endl;
      nFailedFitsAfterCut_++;
    }
  
    return fitres;
}

void PositionFit::finaliseFit(){
  fout_.close();    
  outputFile_->Flush();
  std::cout << " -- Number of invalid fits: " << nInvalidFits_ << std::endl;
  std::cout << " -- Number of fits failed after cutting outliers: " << nFailedFitsAfterCut_ << std::endl;
  std::cout << " -- Number of fits failed for simultaneous fit: reverted back to successful independent fit: " << countFailedSimFits_ << std::endl;

}

unsigned PositionFit::fitEvent(const unsigned ievt, 
			       FitResult & fit,
			       const bool cutOutliers){

  std::vector<unsigned> layerId;
  std::vector<double> posx;
  std::vector<double> posy;
  std::vector<double> posz;
  std::vector<double> posxtruth;
  std::vector<double> posytruth;
  std::vector<double> Ereco;
  layerId.reserve(nLayers_);
  posx.reserve(nLayers_);
  posy.reserve(nLayers_);
  posz.reserve(nLayers_);
  posxtruth.reserve(nLayers_);
  posytruth.reserve(nLayers_);
  Ereco.reserve(nLayers_);
  if (!getPositionFromFile(ievt,
			   layerId,posx,posy,posz,
			   posxtruth,posytruth,
			   Ereco,
			   cutOutliers)){
    nInvalidFits_++;
    return 1;
  }

  const unsigned nL = layerId.size();
  
  if (debug_>0) std::cout << " -- Number of layers for fit: " << nL << std::endl;
  
  //fill some control histograms
  //only once
  if (cutOutliers){
    for (unsigned iL(0); iL<nL;++iL){
      //std::cout << layerId[iL] << " " << posx[iL] << " " << posy[iL] << " " << posz[iL] << std::endl;
      p_recoXvsLayer->Fill(layerId[iL],posx[iL]);
      p_recoYvsLayer->Fill(layerId[iL],posy[iL]);
      p_recoZvsLayer->Fill(layerId[iL],posz[iL]);
      p_truthXvsLayer->Fill(layerId[iL],posxtruth[iL]);
      p_truthYvsLayer->Fill(layerId[iL],posytruth[iL]);
    }
  }

  p_nLayersFit->Fill(nL);

  //if less than 3 valid layers: no point doing a fit !!
  if (nL<3){
    nInvalidFits_++;
    return 1;
  }
  
  
  //Get error matrix removing lines with zero hits
  //TMatrixDSym exx(nL);
  //TMatrixDSym eyy(nL);
  //TMatrixDSym exy(nL);
  //TMatrixDSym eyx(nL);
  TMatrixDSym e(2*nL);
  TVectorD z(2*nL),x(nL),y(nL),xy(2*nL);
  
  for(unsigned i(0);i<nL;++i) {
    z(i)=posz[i];
    z(i+nL)=posz[i];
    if (debug_>0) std::cout << "fit() z(" << i << ") = " << z(i) << std::endl;
    
    for(unsigned j(0);j<nL;++j) {
      //exx(i,j)=matrix_(layerId[i],layerId[j]);
      //eyy(i,j)=matrix_(layerId[i]+nLayers_,layerId[j]+nLayers_);
      //exy(i,j)=matrix_(layerId[i],layerId[j]+nLayers_);
      //eyx(i,j)=matrix_(layerId[i]+nLayers_,layerId[j]);
      e(i,j)=matrix_(layerId[i],layerId[j]);
      e(i+nL,j+nL)=matrix_(layerId[i]+nLayers_,layerId[j]+nLayers_);
      e(i,j+nL)=matrix_(layerId[i],layerId[j]+nLayers_);
      e(i+nL,j)=matrix_(layerId[i]+nLayers_,layerId[j]);

    }
  }

  
  //if (doMatrix_) {
  //exx.Invert();
  // exy.Invert();
  //eyx.Invert();
  //eyy.Invert();
  e.Invert();
  //for(unsigned i(0);i<nL;++i) {
  //for(unsigned j(0);j<nL;++j) {
  // e(i,j) = exx(i,j);
  //e(i+nL,j+nL) = eyy(i,j);
      //e(i+nL,j) = eyx(i,j);
      //e(i,j+nL) = exy(i,j);
  //}
  //}

  //do fit for reco and truth
  //solve:
  // x = tx * (z - z0)
  // y = ty * (z - z0)
  double position0[2][4];
  double positionFF[2][2];
  double position14[2][2];
  double TanAngle[2][2];
 
  for (unsigned rt(0); rt<2;++rt){
    if (debug_) {
      std::cout << "... Processing ";
      if (rt==0) std::cout << " fit to reco position.";
      else std::cout << " fit to truth position.";
      std::cout << std::endl;
    }

    //resolve equation for x and y simultaneously
    if (debug_) {
      std::cout << "... Processing ";
      std::cout << " fit to x-y position.";
      std::cout << std::endl;
    }
    for(unsigned i(0);i<nL;i++) {
      xy(i)= rt==0 ? posx[i] : posxtruth[i];
      x(i) = rt==0 ? posx[i] : posxtruth[i];
      //std::cout << "fit() x(" << i << ") = " << x(i) << std::endl;
    }
    for(unsigned i(0);i<nL;i++) {
      xy(i+nL)= rt==0 ? posy[i] : posytruth[i];
      y(i)=rt==0 ? posy[i] : posytruth[i];
      //std::cout << "fit() x(" << i << ") = " << x(i) << std::endl;
    }

    //first guess
    TMatrixD fitMatrix(4,4);
    TVectorD pars(4);
    for (unsigned ii(0);ii<4;++ii){
      pars(ii) = 0;
      for (unsigned ij(0);ij<4;++ij){
	fitMatrix[ii][ij]=0;
      }
    }
    double chiSq = 0;
    unsigned status = GetIndependentFitResult(e,nL,ievt,
					      z,x,y,xy,
					      fitMatrix,pars,
					      chiSq);
    
    if (status!=0) return status;

    bool doNew = false;
    double ndfSim = 2*nL-3;
    double ndf = 2*nL-4;

    if (!doNew){
      p_chi2[rt]->Fill(chiSq);
      p_chi2overNDF[rt]->Fill(chiSq/ndf);
      if (chiSq/ndf>chi2ndfmax_) return 2;
    }
 

    TVectorD simPars(3);
    TVectorD delta(3);
    //iteration
    unsigned it = 0;
    double chiSqSim = 0;
    double chiSqPrev = chiSq;

    for (unsigned iG(0);iG<3;++iG){//first guess
      simPars(0) = pars(1);
      simPars(1) = pars(3);
      if (iG==0) simPars(2) = -pars(2)/simPars(1);
      else if (iG==1) simPars(2) = -pars(0)/simPars(0);
      else simPars(2) = 0;
      delta(0)=0;
      delta(1)=0;
      delta(2)=0;
      
      if (debug_) {
	std::cout << " ---------------------------------- " << std::endl
		  << " -- Iterations for event " << ievt << std::endl
		  << " ---------------------------------- " << std::endl;
	
	std::cout << " -- First guess = "  << simPars(0) << " " << simPars(1) << " " << simPars(2) << " (" << -pars(0)/simPars(0) << "," << -pars(2)/simPars(1) << ")" << std::endl;
      }
      
      //iteration
      it = 0;
      chiSqSim = 0;
      chiSqPrev = chiSq;
      
      for (; it<100;++it){
	if (debug_) std::cout << " -- iteration " << it << " " << simPars(0) << " " << simPars(1) << " " << simPars(2) << std::endl;
	status = GetSimultaneousFitResult(e,nL,ievt,z,xy,simPars,delta,chiSqSim);
	if (chiSqSim>chiSqPrev && debug_) std::cout << " -- Warning ! Iteration " << it << ", new chi2 = " << chiSqSim << " previous = " << chiSqPrev << std::endl;
	if (debug_) std::cout << " -- deltas = " << delta(0) << " " << delta(1) << " " << delta(2) << std::endl;
	if (chiSqSim<=1.1*chiSq && delta(0) < precision_ && delta(1) < precision_ && delta(2) < precision_) break;
	chiSqPrev = chiSqSim; 
      }
      if (iG<2 && (it==99 || chiSqSim/ndfSim>chi2ndfmax_)){
	if (rt==0 && debug_) std::cout << " -- evt " << ievt << ", failed after iteration " << it << " with chi2=" << chiSqSim << " chi2init = " << chiSq << std::endl;
      }
      else break;
    }
    if (rt==0 && debug_) std::cout << " -- evt " << ievt << ", found result after iteration " << it << " with chi2=" << chiSqSim << " chi2init = " << chiSq << std::endl;
    
    p_iterations->Fill(it);
    if (doNew && status!=0) return status;
    
    if (doNew) {
      //if not found, revert back to independent fit result...
      if (it==99 || chiSqSim/ndfSim>chi2ndfmax_){
	if (chiSq/ndf>chi2ndfmax_) return 2;
	else {
	  chiSqSim = chiSq;
	  ndfSim = ndf;
	  doNew = false;
	  if (rt==0) countFailedSimFits_++;
	}
      }
      p_chi2[rt]->Fill(chiSqSim);
      p_chi2overNDF[rt]->Fill(chiSqSim/ndfSim);
    }
    
    if (doNew){
      position0[rt][0] = -1.*simPars(0)*simPars(2);
      position0[rt][1] = -1.*simPars(1)*simPars(2);
      position0[rt][2] = simPars(2);
      position0[rt][3] = simPars(2);
      positionFF[rt][0] = simPars(0)*(posz[0]-simPars(2));
      positionFF[rt][1] = simPars(1)*(posz[0]-simPars(2));
      position14[rt][0] = simPars(0)*(posz[14]-simPars(2));
      position14[rt][1] = simPars(1)*(posz[14]-simPars(2));
      TanAngle[rt][0] = simPars(0);
      TanAngle[rt][1] = simPars(1);
    }
    else {
      position0[rt][0] = pars(0);
      position0[rt][1] = pars(2);
      position0[rt][2] = -1*pars(0)/pars(1);
      position0[rt][3] = -1*pars(2)/pars(3);
      positionFF[rt][0] = pars(0)+pars(1)*posz[0];
      positionFF[rt][1] = pars(2)+pars(3)*posz[0];
      position14[rt][0] = pars(0)+pars(1)*posz[14];
      position14[rt][1] = pars(2)+pars(3)*posz[14];
      TanAngle[rt][0] = pars(1);
      TanAngle[rt][1] = pars(3);
    }

    p_impactX0[rt]->Fill(position0[rt][0]);
    p_impactXFF[rt]->Fill(positionFF[rt][0]);
    p_impactX14[rt]->Fill(position14[rt][0]);
    p_tanAngleX[rt]->Fill(TanAngle[rt][0]);
    p_impactY0[rt]->Fill(position0[rt][1]);
    p_impactYFF[rt]->Fill(positionFF[rt][1]);
    p_impactY14[rt]->Fill(position14[rt][1]);
    p_tanAngleY[rt]->Fill(TanAngle[rt][1]);
    p_impactZx[rt]->Fill(position0[rt][2]);
    p_impactZy[rt]->Fill(position0[rt][3]);

    if (rt==0) {
      p_positionReso->Fill(sqrt(fabs(fitMatrix[0][0])));
      p_angularReso->Fill(sqrt(fabs(fitMatrix[1][1])));
      
      fout_ << ievt << " " 
	    << position0[0][0] << " " 
	    << sqrt(fabs(fitMatrix[0][0])) << " " 
	    << TanAngle[0][0] << " " 
	    << sqrt(fabs(fitMatrix[1][1])) << " "
	    << position0[0][1] << " " 
	    << sqrt(fabs(fitMatrix[2][2])) << " "
	    << TanAngle[0][1] << " "
	    << sqrt(fabs(fitMatrix[3][3]));
      //<< std::endl;

      for (unsigned iL(0); iL<nL;++iL){
	double x = position0[0][0]+TanAngle[0][0]*posz[iL];
	double y = position0[0][1]+TanAngle[0][1]*posz[iL];
	p_fitXvsLayer->Fill(layerId[iL],x);
	p_fitYvsLayer->Fill(layerId[iL],y);
	//eventPos[layerId[iL]] = ROOT::Math::XYZVector(x,y,posz[iL]);
      }
      fit.found = true;
      fit.pos_x = position0[0][0];
      fit.tanangle_x = TanAngle[0][0];
      fit.pos_y = position0[0][1];
      fit.tanangle_y = TanAngle[0][1];
    }//reco
    else {
      fout_ << " " << position0[1][0]
	    << " " << TanAngle[1][0]
	    << " " << position0[1][1]
	    << " " << TanAngle[1][1]
	    << std::endl;
    }//truth

  }//reco or truth

  recoDir_ = Direction(TanAngle[0][0],TanAngle[0][1]);

  p_impactX0_residual->Fill(position0[0][0]-position0[1][0]);
  p_impactXFF_residual->Fill(positionFF[0][0]-positionFF[1][0]);
  p_impactX14_residual->Fill(position14[0][0]-position14[1][0]);
  p_tanAngleX_residual->Fill(TanAngle[0][0]-TanAngle[1][0]);
  p_angleX_residual->Fill(atan(TanAngle[0][0])-atan(TanAngle[1][0]));
  p_impactY0_residual->Fill(position0[0][1]-position0[1][1]);
  p_impactYFF_residual->Fill(positionFF[0][1]-positionFF[1][1]);
  p_impactY14_residual->Fill(position14[0][1]-position14[1][1]);
  p_impactZx_residual->Fill(position0[0][2]-position0[1][2]);
  p_impactZy_residual->Fill(position0[0][3]-position0[1][3]);
  p_tanAngleY_residual->Fill(TanAngle[0][1]-TanAngle[1][1]);
  p_angleY_residual->Fill(atan(TanAngle[0][1])-atan(TanAngle[1][1]));
      
  Direction truthDir = Direction(TanAngle[1][0],TanAngle[1][1]);

  p_eta_reco->Fill(recoDir_.eta());
  p_phi_reco->Fill(recoDir_.phi());
  p_eta_truth->Fill(truthDir.eta());
  p_phi_truth->Fill(truthDir.phi());
  p_eta_residual->Fill(recoDir_.eta()-truthDir.eta());
  p_phi_residual->Fill(recoDir_.phi()-truthDir.phi());

  //std::cout << " -- Size of eventPos=" << eventPos.size() << std::endl;
  return 0;
}

unsigned PositionFit::GetIndependentFitResult(const TMatrixDSym & e,
					      const unsigned nL,
					      //const unsigned rt,
					      const unsigned ievt,
					      const TVectorD & z,
					      const TVectorD & x,
					      const TVectorD & y,
					      const TVectorD & xy,
					      TMatrixD & w,
					      TVectorD & p,
					      double & chiSq
					      ){

  TVectorD v(4);
  //TVectorD vx(2);
  //TVectorD px(2);
  //TVectorD vy(2);
  //TVectorD py(2);
  TVectorD ux(2*nL);
  TVectorD uy(2*nL);
  TVectorD zx(2*nL);
  TVectorD zy(2*nL);
  //TVectorD ux(nL);
  //TMatrixD wx(2,2);
  // TMatrixD wy(2,2);
  //TMatrixDSym ex(nL);
  //TMatrixDSym ey(nL);

  for(unsigned i(0);i<nL;++i) {
    ux(i)=1.0;
    uy(i)=0.0;
    uy(i+nL)=1.0;
    ux(i+nL)=0.0;
    zx(i)=z(i);
    zx(i+nL)=0;
    zy(i+nL)=z(i);
    zy(i)=0;
  }

  w(0,0)=ux*(e*ux);
  w(0,1)=ux*(e*zx);
  w(0,2)=ux*(e*uy);
  w(0,3)=ux*(e*zy);

  w(1,0)=zx*(e*ux);
  w(1,1)=zx*(e*zx);
  w(1,2)=zx*(e*uy);
  w(1,3)=zx*(e*zy);

  w(2,0)=uy*(e*ux);
  w(2,1)=uy*(e*zx);
  w(2,2)=uy*(e*uy);
  w(2,3)=uy*(e*zy);

  w(3,0)=zy*(e*ux);
  w(3,1)=zy*(e*zx);
  w(3,2)=zy*(e*uy);
  w(3,3)=zy*(e*zy);

  w.Invert();

  v(0)=ux*(e*xy);
  v(1)=zx*(e*xy);
  v(2)=uy*(e*xy);
  v(3)=zy*(e*xy);

  p=w*v;


  /*for(unsigned i(0);i<nL;++i) {
    ux(i)=1.0;
    zx(i)=z(i);
    for(unsigned j(0);j<nL;++j) {
      ex(i,j)=e(i,j);
      ey(i,j)=e(i+nL,j+nL);
    }
    }

  wx(0,0)=ux*(ex*ux);
  wx(0,1)=ux*(ex*zx);
  wx(1,0)=zx*(ex*ux);
  wx(1,1)=zx*(ex*zx);

  wy(0,0)=ux*(ey*ux);
  wy(0,1)=ux*(ey*zx);
  wy(1,0)=zx*(ey*ux);
  wy(1,1)=zx*(ey*zx);
 
  vx(0)=ux*(ex*x);
  vx(1)=zx*(ex*x);
  vy(0)=ux*(ey*y);
  vy(1)=zx*(ey*y);
  
  wx.Invert();
  wy.Invert();
  for(unsigned i(0);i<2;++i) {
    for(unsigned j(0);j<2;++j) {
      w(i,j)=wx(i,j);
      w(i+2,j+2)=wy(i,j);
    }
  }
  px=wx*vx;
  py=wy*vy;
  p(0)=px(0);
  p(1)=px(1);
  p(2)=py(0);
  p(3)=py(1);
  */

  if (debug_) {
    std::cout << "fit() w(0,0) = " << w(0,0) << std::endl;
    std::cout << "fit() w(0,1) = " << w(0,1) << std::endl;
    std::cout << "fit() w(1,0) = " << w(1,0) << std::endl;
    std::cout << "fit() w(1,1) = " << w(1,1) << std::endl;	
    std::cout << "fit() w(2,2) = " << w(2,2) << std::endl;
    std::cout << "fit() w(2,3) = " << w(2,3) << std::endl;
    std::cout << "fit() w(3,2) = " << w(3,2) << std::endl;
    std::cout << "fit() w(3,3) = " << w(3,3) << std::endl;	
    std::cout << "fit() p(0) = " << p(0) << std::endl;
    std::cout << "fit() p(1) = " << p(1) << std::endl;
    std::cout << "fit() p(2) = " << p(2) << std::endl;
    std::cout << "fit() p(3) = " << p(3) << std::endl;
  }
  
  TVectorD dp(2*nL);
  //TVectorD dpx(nL);
  //TVectorD dpy(nL);
  for(unsigned i(0);i<nL;i++) {
    dp(i)=xy(i)-p(0)-p(1)*z(i);
    dp(i+nL)=xy(i+nL)-p(2)-p(3)*z(i);
    //dpx(i)=x(i)-px(0)-px(1)*zx(i);
    //dpy(i)=y(i)-py(0)-py(1)*zx(i);
  }
  
  //double chiSq=dpx*(ex*dpx)+dpy*(ey*dpy);
  chiSq=dp*(e*dp);

  //chi2 test
  //number of points: x and y per layer minus number of parameters
  double ndf = 2*nL-4;

  if (chiSq/ndf>chi2ndfmax_) {
      std::cout << " ---- Independent Fit failed for event " << ievt << std::endl;
      std::cout << "Chi2/ndf = " << chiSq << "/" << ndf << "=" << chiSq/ndf << std::endl;
      //std::cout << "fitw(0,0) = " << fitMatrix[0][0] << std::endl;
      //std::cout << "fitw(1,1) = " << fitMatrix[1][1] << std::endl;
      //std::cout << "fitw(2,2) = " << fitMatrix[2][2] << std::endl;
      //std::cout << "fitw(3,3) = " << fitMatrix[3][3] << std::endl;	
      std::cout << "position = " << p(0) << " " << p(2) << std::endl;
      std::cout << "tanAngle = " << p(1) << " " << p(3) << std::endl;
      //return 2;
  }
  

  return 0;

}

unsigned PositionFit::GetSimultaneousFitResult(const TMatrixDSym & e,
					       const unsigned nL,
					       //const unsigned rt,
					       const unsigned ievt,
					       const TVectorD & z,
					       const TVectorD & xy,
					       TVectorD & simPars,
					       TVectorD & delta,
					       double & chiSq
					       ){
  //first derivative using first guess
  TVectorD yminusf(2*nL);
  TVectorD f(2*nL),dfdp0(2*nL),dfdp1(2*nL),dfdp2(2*nL);
  std::vector<TVectorD> dfdp;
  TVectorD d2fdp0dp2(2*nL),d2fdp1dp2(2*nL);
  std::vector<std::vector<TVectorD> > d2fdp;
  for(unsigned i(0);i<2*nL;++i) {
    if (i<nL) {
      f(i) = simPars(0)*(z(i)-simPars(2));
      dfdp0(i) = (z(i)-simPars(2));
      dfdp1(i) = 0;
      dfdp2(i) = -1.*simPars(0);
      d2fdp0dp2(i) = -1;
      d2fdp1dp2(i) = 0;
    }
    else {
      f(i) = simPars(1)*(z(i)-simPars(2));
      dfdp0(i) = 0;
      dfdp1(i) = (z(i)-simPars(2));
      dfdp2(i) = -1.*simPars(1);
      d2fdp0dp2(i) = 0;
      d2fdp1dp2(i) = -1;
    }
    yminusf(i) = xy(i)-f(i);
  }
  dfdp.push_back(dfdp0);
  dfdp.push_back(dfdp1);
  dfdp.push_back(dfdp2);
  std::vector<TVectorD> empty;
  TVectorD emptyvec(2*nL);
  for(unsigned i(0);i<2*nL;++i) emptyvec(i)=0;
  empty.resize(3,emptyvec);
  d2fdp.resize(3,empty);
  d2fdp[0][2] = d2fdp0dp2;
  d2fdp[1][2] = d2fdp1dp2;
  d2fdp[2][0] = d2fdp0dp2;
  d2fdp[2][1] = d2fdp1dp2;
  
  TMatrixD G(3,3);
  TVectorD g0(3);
  for(unsigned i(0);i<3;++i) {
    g0(i) = -2.*yminusf*(e*dfdp[i]);
    for(unsigned j(0);j<3;++j) {
      G(i,j) = -2.*(-1.*dfdp[i]*(e*dfdp[j]) + yminusf*(e*d2fdp[i][j]));
    }
  }
  G.Invert();
  
  delta = -1.*G*g0;
  
  simPars += delta;
  
  for(unsigned i(0);i<2*nL;++i) {
    if (i<nL)
      f(i) = simPars(0)*(z(i)-simPars(2));
    else
      f(i) = simPars(1)*(z(i)-simPars(2));
    yminusf(i) = xy(i)-f(i);
  }
  
  chiSq=yminusf*(e*yminusf);
  //double ndf = 2*nL-3;
  
  /*if (chiSq/ndf>chi2ndfmax_) {
    std::cout << " ---- Fit failed for event " << ievt << std::endl;
    std::cout << "Chi2/ndf = " << chiSq << "/" << ndf << "=" << chiSq/ndf << std::endl;
    std::cout << "position = " << simPars(2) << std::endl;
    std::cout << "tanAngle = " << simPars(0) << " " << simPars(1) << std::endl;
    return 2;
    }*/
  
  return 0;
  
}


bool PositionFit::getPositionFromFile(const unsigned ievt,
				      std::vector<unsigned> & layerId,
				      std::vector<double> & posx,
				      std::vector<double> & posy,
				      std::vector<double> & posz,
				      std::vector<double> & posxtruth,
				      std::vector<double> & posytruth,
				      std::vector<double> & E,
				      bool cutOutliers,
				      bool print){


  std::ifstream fin;
  std::ostringstream finname;
  finname << outFolder_ << "/initialPos_evt" << ievt << ".dat";
  fin.open(finname.str());
  if (!fin.is_open()){
    if (print) std::cout << " Cannot open input file " << finname.str() << "!" << std::endl;
      return false;
  }
  
  while (!fin.eof()){
    unsigned l=nLayers_;
    double xr=0,yr=0,xt=0,yt=0,e=0;
    fin>>l>>xr>>yr>>xt>>yt;
    if (!doMatrix_) fin>>e;
    if (l<nLayers_){
      //bool l7to22 = true;//l>6 && l<23;
      bool pass = fabs(xr-xt)<residualMax_ && fabs(yr-yt)<residualMax_;
      //unsigned posmm = static_cast<unsigned>(fabs(yt)+5);
      //bool isEdge = true;//posmm%10 <= 2 || posmm%10 >= 8;
      if (!cutOutliers || (cutOutliers && pass)){
	layerId.push_back(l);
	posx.push_back(xr);
	posy.push_back(yr);
	posz.push_back(avgZ_[l]);
	posxtruth.push_back(xt);
	posytruth.push_back(yt);
      }
      //use all for energy estimate
      if (!doMatrix_) E.push_back(e);
    }
  }
  
  fin.close();
  /*
  //@TODO to use something else than truth info :/
  if (cutOutliers){
    //
    for (unsigned i(0);i<layerId.size();++i){
      double xt=?;
      double yt=?;
      double xr=posx[i];
      double yr=posy[i];
      bool pass=fabs(xr-xt)<residualMax_ && fabs(yr-yt)<residualMax_;
      std::cout << i << " l=" << layerId[i] << " pass=" << pass << " xr=" << xr << " yr=" << yr << std::endl;
      if (!pass) {
	std::cout << " ---- erase layer " << *(layerId.begin()+i) << std::endl;
	layerId.erase(layerId.begin()+i);
	posx.erase(posx.begin()+i);
	posy.erase(posy.begin()+i);
	posxtruth.erase(posxtruth.begin()+i);
	posytruth.erase(posytruth.begin()+i);
	i--;
      }
    }
  }
  */

  return true;
}




