#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#include <boost/algorithm/string.hpp>
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "boost/function.hpp"

// helpful tools
#include "KDTreeLinkerAlgoT.h"
#include <unordered_map>
#include <unordered_set>

#include "TFile.h"

#include "HGCSSEvent.hh"


typedef KDTreeLinkerAlgo<unsigned,3> KDTree;
typedef KDTreeNodeInfoT<unsigned,3> KDNode;


namespace {

  std::pair<float,float> minmax(const float a, const float b) {
    return ( b < a ? 
	     std::pair<float,float>(b, a) : 
	     std::pair<float,float>(a, b)   );
  }
  
  template<typename T>
  KDTreeCube fill_and_bound_kd_tree(const std::vector<T>& points,
				    const std::vector<bool>& usable,
				    std::vector<KDTreeNodeInfoT<unsigned,3> >& nodes) {
    std::array<float,3> minpos{ {0.0f,0.0f,0.0f} }, maxpos{ {0.0f,0.0f,0.0f} };
    for( unsigned i = 0 ; i < points.size(); ++i ) {
      if( !usable[i] ) continue;
      const auto& pos = points[i].position();
      nodes.emplace_back(i, (float)pos.X(), (float)pos.Y(), (float)pos.Z());

      //std::cout << " ele " << i << " position: " << pos.X() << " " << pos.Y() << " " << pos.Z() << " nodes size:" << nodes.size() << std::endl;

      if( i == 0 ) {
	minpos[0] = pos.X(); minpos[1] = pos.Y(); minpos[2] = pos.Z();
	maxpos[0] = pos.X(); maxpos[1] = pos.Y(); maxpos[2] = pos.Z();
      } else {
	minpos[0] = std::min((float)pos.X(),minpos[0]);
	minpos[1] = std::min((float)pos.Y(),minpos[1]);
	minpos[2] = std::min((float)pos.Z(),minpos[2]);
	maxpos[0] = std::max((float)pos.X(),maxpos[0]);
	maxpos[1] = std::max((float)pos.Y(),maxpos[1]);
	maxpos[2] = std::max((float)pos.Z(),maxpos[2]);
      }
    }
    
    //std::cout << " KDTree filled ! -- nElements: " << nodes.size() << std::endl;
    //std::cout << " -- min and max X = " << minpos[0] << " " << maxpos[0] << std::endl
    //<< " -- min and max Y = " << minpos[1] << " " << maxpos[1] << std::endl
    //<< " -- min and max Z = " << minpos[2] << " " << maxpos[2] << std::endl;

    return KDTreeCube(minpos[0],maxpos[0],
		      minpos[1],maxpos[1],
		      minpos[2],maxpos[2]);
  }
}//namespace


bool testInputFile(std::string input, TFile* & file){
  file = TFile::Open(input.c_str());
  
  if (!file) {
    std::cout << " -- Error, input file " << input.c_str() << " cannot be opened. Skipping..." << std::endl;
    return false;
  }
  else std::cout << " -- input file " << file->GetName() << " successfully opened." << std::endl;
  return true;
};

double zPos(const int l){
  if (l<=0) return 3173.9;
  if (l==1) return 3183.65;
  if (l==2) return 3190.95;
  if (l==3) return 3200.4;
  if (l==4) return 3207.7;
  if (l==5) return 3217.15;
  if (l==6) return 3224.45;
  if (l==7) return 3233.9;
  if (l==8) return 3241.2;
  if (l==9) return 3250.65;
  if (l==10) return 3257.95;
  if (l==11) return 3268.45;
  if (l==12) return 3276.85;
  if (l==13) return 3287.35;
  if (l==14) return 3295.75;
  if (l==15) return 3306.25;
  if (l==16) return 3314.65;
  if (l==17) return 3325.15;
  if (l==18) return 3333.55;
  if (l==19) return 3344.05;
  if (l==20) return 3352.45;
  if (l==21) return 3364.35;
  if (l==22) return 3375.05;
  if (l==23) return 3386.95;
  if (l==24) return 3397.65;
  if (l==25) return 3409.55;
  if (l==26) return 3420.25;
  if (l==27) return 3432.15;
  if (l==28) return 3442.85;
  if (l>=29) return 3454.75;
  return 0;
};

double calibratedE(const double Etot, const double eta){
  //calibration for signal region 2: 3*3 cm^2
  double pars[3] = {77,3.4,-0.50};
  double paro[3] = {-11.6,-7.7,-8.8};
  //double paro[3] = {-5.3,-12.8,-6.9};  
  double offset = paro[0] + paro[1]*eta + paro[2]*eta*eta;
  double slope = pars[0] + pars[1]*eta + pars[2]*eta*eta;
  return (Etot-offset)/slope;
};

double DeltaPhi(const double & phi1, const double & phi2){
  double dphi = phi1 - phi2;
  if (dphi< (-1.*TMath::Pi())) dphi += 2*TMath::Pi();
  if (dphi>TMath::Pi()) dphi -= 2*TMath::Pi();
  return dphi;
};

double E(const unsigned pT, const unsigned eta){
  return pT*cosh(eta/10.);
};

double pT(const unsigned E, const unsigned eta){
  return E/cosh(eta/10.);
};

double absWeight(const unsigned layer, const bool dedx=true){
  if (dedx==false){
    if (layer == 0) return 0.0378011;//6.285/95.4=0.06588;
    if (layer == 1) return 1;//95.4/95.4=1
    if (layer == 2) return 0.646989;//88.16/95.4=0.92
    if (layer == 3) return 0.617619;//51.245/95.4=0.537
    if (layer == 4) return 0.646989;
    if (layer == 5) return 0.617619;
    if (layer == 6) return 0.646989;
    if (layer == 7) return 0.617619;
    if (layer == 8) return 0.646989;
    if (layer == 9) return 0.617619;
    if (layer == 10) return 0.646989;
    if (layer == 11) return 0.942829;//74.45/95.4=0.78
    if (layer == 12) return 0.859702;//102.174/95.4=1.071
    if (layer == 13) return 0.942829;
    if (layer == 14) return 0.859702;
    if (layer == 15) return 0.942829;
    if (layer == 16) return 0.859702;
    if (layer == 17) return 0.942829;
    if (layer == 18) return 0.859702;
    if (layer == 19) return 0.942829;
    if (layer == 20) return 0.859702;
    if (layer == 21) return 1.37644;//105.39/95.4=1.1047
    if (layer == 22) return 1.30447;//131.476/95.4=1.378
    if (layer == 23) return 1.37644;
    if (layer == 24) return 1.30447;
    if (layer == 25) return 1.37644;
    if (layer == 26) return 1.30447;
    if (layer == 27) return 1.37644;
    if (layer == 28) return 1.30447;
    if (layer == 29) return 1.37644;//1.79662;//
  }
  else {
    if (layer == 0) return 0.248;
    if (layer == 1) return 1;//95.4/95.4=1
    if (layer == 2) return 0.92;//88.16/95.4=0.92
    if (layer == 3) return 0.537;//51.245/95.4=0.537
    if (layer == 4) return 0.92;
    if (layer == 5) return 0.537;
    if (layer == 6) return 0.92;
    if (layer == 7) return 0.537;
    if (layer == 8) return 0.92;
    if (layer == 9) return 0.537;
    if (layer == 10) return 0.92;
    if (layer == 11) return 0.78;//74.45/95.4=0.78
    if (layer == 12) return 1.071;//102.174/95.4=1.071
    if (layer == 13) return 0.78;
    if (layer == 14) return 1.071;
    if (layer == 15) return 0.78;
    if (layer == 16) return 1.071;
    if (layer == 17) return 0.78;
    if (layer == 18) return 1.071;
    if (layer == 19) return 0.78;
    if (layer == 20) return 1.071;
    if (layer == 21) return 1.1047;//105.39/95.4=1.1047
    if (layer == 22) return 1.378;//131.476/95.4=1.378
    if (layer == 23) return 1.1047;
    if (layer == 24) return 1.378;
    if (layer == 25) return 1.1047;
    if (layer == 26) return 1.378;
    if (layer == 27) return 1.1047;
    if (layer == 28) return 1.378;
    if (layer == 29) return 1.1047;
  }
  return 1;
};

double getEtaPhys(const double & posx,
		  const double & posy,
		  const double & posz,
		  const HGCSSEvent * event){

    double posxcor = posx-event->vtx_x();
    double posycor = posy-event->vtx_y();
    double poszcor = posz-event->vtx_z();
    double r3d = sqrt(pow(posxcor,2)+pow(posycor,2)+pow(poszcor,2));
    double thetaPhys = acos(fabs(poszcor)/r3d);
    double etaPhys = (poszcor)>0? -1.*log(tan(thetaPhys/2.)) : log(tan(thetaPhys/2.));
    return etaPhys;
};
