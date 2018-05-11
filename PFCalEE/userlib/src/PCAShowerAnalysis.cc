//#include "RecoEgamma/Examples/interface/PCAShowerAnalysis.h"
#include "PCAShowerAnalysis.h"

#include "TMatrixDfwd.h"
#include "TVectorD.h"

PCAShowerAnalysis::PCAShowerAnalysis ( bool segmented, 
				       bool logweighting, 
				       bool debug) : 
  principal_(0), 
  logweighting_(logweighting), 
  segmented_(segmented), 
  alreadyfilled_(false), 
  debug_(debug)
{

  principal_ = new TPrincipal(3,"D"); 

  // minimal rechit value
  mip_ = 0.000055;//40;
  entryz_ = 320.38;
  debug_=false;
}

PCAShowerAnalysis::~PCAShowerAnalysis ()
{
  delete principal_;
}

void PCAShowerAnalysis::showerParameters(const HGCSSCluster & clus)
{

  if (!alreadyfilled_) {
    
    double variables[3] = {0.,0.,0.};
    const std::map<HGCSSRecoHit*,double> & lmap = clus.recHitFractions();
    std::map<HGCSSRecoHit*,double>::const_iterator iter = lmap.begin();
    
    if (debug_) std::cout << " -- Number of rechits in cluster: " << clus.nRecHits() << " " << lmap.size() << std::endl;
    unsigned counter = 0;
    for (;iter!=lmap.end();++iter){
      HGCSSRecoHit* myhit = iter->first;
      if (!myhit) {
	if (debug_) std::cout << " Hit " << myhit << " not found..." << std::endl;
	continue;
      }
      ROOT::Math::XYZPoint cellPos(myhit->position());
      double en = myhit->energy();
      variables[0] = cellPos.x(); 
      variables[1] = cellPos.y(); 
      variables[2] = cellPos.z();
      if (!segmented_) variables[2] = entryz_;
      if (!logweighting_) {
	// energy weighting
	for (int i=0; i<int(en); i++) principal_->AddRow(variables); 
      } else {
	// a log-weighting, energy not in fraction of total
	double w0 = -log(20.); // threshold, could use here JB's thresholds
	double scale = 250.; // to scale the weight so to get ~same nbr of points as for E-weight 
	                     //  for the highest hit of ~0.1 GeV
	int nhit = int(scale*(w0+log(en)));
	if (nhit<=0) nhit=1;
	for (int i=0; i<nhit; i++) principal_->AddRow(variables);	
      }
      counter++;
      //if (debug_) std::cout << " - hit added n=" << counter << std::endl;
    }
    if (counter!=lmap.size()) std::cout << " -- Warning, not all hits found for making principals ! Found " << counter << " out of " << lmap.size() << std::endl;
  }


  if (debug_) std::cout << " Making principals " << std::endl;

  alreadyfilled_ = true;
  principal_->MakePrincipals();
  
  if (debug_) std::cout << "*** Principal component analysis (standalone) ****" << std::endl;
  if (debug_) std::cout << "shower average (x,y,z) = " << "(" << (*principal_->GetMeanValues())[0] << ", " <<
		(*principal_->GetMeanValues())[1] << ", " << (*principal_->GetMeanValues())[2] << ")" << std::endl;

  TMatrixD matrix = *principal_->GetEigenVectors();
  if (debug_) std::cout << "shower main axis (x,y,z) = " << "(" << matrix(0,0) << ", " <<
		matrix(1,0) << ", " << matrix(2,0) << ")" << std::endl;

  TVectorD eigenvalues = *principal_->GetEigenValues();
  if (debug_) std::cout << "shower eigen values = " 
			<< "(" << eigenvalues(0) << ", " 
			<< eigenvalues(1) << ", " 
			<< eigenvalues(2) << ")" 
			<< std::endl;

  TVectorD sigmas = *principal_->GetSigmas();
  if (debug_) std::cout << "shower sigmas = " << "(" << sigmas(0) << ", " <<
		sigmas(1) << ", " << sigmas(2) << ")" << std::endl;
  


  showerBarycenter = ROOT::Math::XYZPoint((*principal_->GetMeanValues())[0],(*principal_->GetMeanValues())[1],(*principal_->GetMeanValues())[2]);	 
  showerAxis = ROOT::Math::XYZVector(matrix(0,0),matrix(1,0),matrix(2,0));
   
  showerEigenValues = ROOT::Math::XYZVector(eigenvalues(0), eigenvalues(1), eigenvalues(2));

  showerSigmas = ROOT::Math::XYZVector(sigmas(0), sigmas(1), sigmas(2));


  // resolve direction ambiguity
  if (showerAxis.z()*showerBarycenter.z()<0) {
    showerAxis = ROOT::Math::XYZVector(-matrix(0,0),-matrix(1,0),-matrix(2,0));
    if (debug_) std::cout << "PCA shower dir reverted " << showerAxis << "eta " << showerAxis.eta() << " phi " << showerAxis.phi() << std::endl;
  }
	 
  return;

}

  
