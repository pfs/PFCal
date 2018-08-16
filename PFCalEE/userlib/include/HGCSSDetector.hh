#ifndef HGCSSDetector_h
#define HGCSSDetector_h

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "TH2D.h"

enum DetectorEnum {
  FECAL,
  MECAL,
  BECAL,
  FHCAL,
  BHCAL1,
  BHCAL2
};

class HGCSSSubDetector {

public:
  HGCSSSubDetector():
    type(FECAL),
    name(""),
    layerIdMin(0),
    layerIdMax(0),
    mipWeight(1),
    absWeight(1),
    gevWeight(1),
    gevOffset(0),
    isSi(false),
    isScint(false),
    radiusLim(0)
  {};
  ~HGCSSSubDetector(){};

  DetectorEnum type;
  std::string name;
  unsigned layerIdMin;
  unsigned layerIdMax;
  double mipWeight;
  double absWeight;
  double gevWeight;
  double gevOffset;
  bool isSi;
  bool isScint;
  double radiusLim;

  inline unsigned nLayers() const{
    return (layerIdMax-layerIdMin);
  };

private:

};

class HGCSSDetector {

public:
  friend HGCSSDetector & theDetector();

  inline void initialiseIndices(const unsigned versionNumber){
    
    indices_.clear();
    indices_.resize(7,0);
    //fill layer indices
    if (versionNumber==22){
      indices_[4] = 0;
      indices_[5] = 10;
      indices_[6] = 10;
    }
    else if (versionNumber==28 || versionNumber==32) {
      indices_[4] = 0;
      indices_[5] = 12;
      indices_[6] = 12;
    }
    else if (versionNumber==23) {
      indices_[3] = 0;
      indices_[4] = 38;
      indices_[5] = 47;
      indices_[6] = 54;
    }
    else if (versionNumber==21) {
      indices_[3] = 0;
      indices_[4] = 24;
      indices_[5] = 34;
      indices_[6] = 34;
    }
    else if (versionNumber==27 || versionNumber==31) {
      indices_[3] = 0;
      indices_[4] = 12;
      indices_[5] = 24;
      indices_[6] = 24;
    }
    else if (versionNumber==38) {
      indices_[3] = 0;
      indices_[4] = 11;
      indices_[5] = 23;
      indices_[6] = 23;
    }
    else if (versionNumber==39) {
      indices_[3] = 0;
      indices_[4] = 9;
      indices_[5] = 21;
      indices_[6] = 21;
    }
    else if (versionNumber < 20){
      indices_[0] = 0;
      indices_[1] = versionNumber==8?11:10;
      indices_[2] = versionNumber==8?21:20;
      indices_[3] = versionNumber==8?31:30;
      indices_[4] = indices_[3];
      indices_[5] = indices_[3];
      indices_[6] = indices_[3];
    }
    else if (versionNumber == 30 || versionNumber == 60 || (versionNumber >= 100 && versionNumber < 104)){
      indices_[0] = 0;
      indices_[1] = 10;
      indices_[2] = 20;
      indices_[3] = 28;
      indices_[4] = indices_[3];
      indices_[5] = indices_[3];
      indices_[6] = indices_[3];

      sensitiveZ_.resize(indices_[3],0);
      sensitiveZ_[0] = 3198;
      sensitiveZ_[1] = 3207.1;
      sensitiveZ_[2] = 3222.4;
      sensitiveZ_[3] = 3231.5;
      sensitiveZ_[4] = 3246.8;
      sensitiveZ_[5] = 3255.9;
      sensitiveZ_[6] = 3271.2;
      sensitiveZ_[7] = 3280.3;
      sensitiveZ_[8] = 3295.6;
      sensitiveZ_[9] = 3304.7;
      sensitiveZ_[10] = 3320;
      sensitiveZ_[11] = 3329.1;
      sensitiveZ_[12] = 3344.4;
      sensitiveZ_[13] = 3353.5;
      sensitiveZ_[14] = 3368.8;
      sensitiveZ_[15] = 3377.9;
      sensitiveZ_[16] = 3393.2;
      sensitiveZ_[17] = 3402.3;
      sensitiveZ_[18] = 3417.6;
      sensitiveZ_[19] = 3426.7;
      sensitiveZ_[20] = 3442;
      sensitiveZ_[21] = 3451.1;
      sensitiveZ_[22] = 3466.4;
      sensitiveZ_[23] = 3475.5;
      sensitiveZ_[24] = 3490.8;
      sensitiveZ_[25] = 3499.9;
      sensitiveZ_[26] = 3515.2;
      sensitiveZ_[27] = 3524.3;

    }
    else if (versionNumber == 33){
      indices_[0] = 0;
      indices_[1] = 10;
      indices_[2] = 20;
      indices_[3] = 28;
      indices_[4] = 40;
      indices_[5] = 52;
      indices_[6] = 52;
    }
    else if (versionNumber == 34){
      indices_[0] = 0;
      indices_[1] = 8;
      indices_[2] = 16;
      indices_[3] = 24;
      indices_[4] = indices_[3];
      indices_[5] = indices_[3];
      indices_[6] = indices_[3];
    }
    else if (versionNumber == 35){
      indices_[0] = 0;
      indices_[1] = 6;
      indices_[2] = 12;
      indices_[3] = 18;
      indices_[4] = indices_[3];
      indices_[5] = indices_[3];
      indices_[6] = indices_[3];
    }
    else if (versionNumber == 36){
      indices_[0] = 0;
      indices_[1] = 8;
      indices_[2] = 16;
      indices_[3] = 24;
      indices_[4] = 35;
      indices_[5] = 47;
      indices_[6] = 47;
    }
    else if (versionNumber == 37){
      indices_[0] = 0;
      indices_[1] = 6;
      indices_[2] = 12;
      indices_[3] = 18;
      indices_[4] = 27;
      indices_[5] = 39;
      indices_[6] = 39;
    }
    else if (versionNumber == 61){
      indices_[3] = 0;
      indices_[4] = 25;
      indices_[5] = 41;
      indices_[6] = 41;

      sensitiveZ_.resize(indices_[6],0);
      sensitiveZ_[0] = 3037;
      sensitiveZ_[1] = 3086;
      sensitiveZ_[2] = 3135;
      sensitiveZ_[3] = 3184;
      sensitiveZ_[4] = 3233;
      sensitiveZ_[5] = 3282;
      sensitiveZ_[6] = 3331;
      sensitiveZ_[7] = 3380;
      sensitiveZ_[8] = 3429;
      sensitiveZ_[9] = 3479.9;
      sensitiveZ_[10] = 3530.8;
      sensitiveZ_[11] = 3581.7;
      sensitiveZ_[12] = 3665.6;
      sensitiveZ_[13] = 3749.5;
      sensitiveZ_[14] = 3833.4;
      sensitiveZ_[15] = 3917.3;
      sensitiveZ_[16] = 4001.2;
      sensitiveZ_[17] = 4085.1;
      sensitiveZ_[18] = 4169;
      sensitiveZ_[19] = 4252.9;
      sensitiveZ_[20] = 4336.8;
      sensitiveZ_[21] = 4420.7;
      sensitiveZ_[22] = 4504.6;
      sensitiveZ_[23] = 4588.5;
      sensitiveZ_[24] = 0;
      sensitiveZ_[25] = 3430.8;
      sensitiveZ_[26] = 3481.7;
      sensitiveZ_[27] = 3532.6;
      sensitiveZ_[28] = 3583.5;
      sensitiveZ_[29] = 3667.4;
      sensitiveZ_[30] = 3751.3;
      sensitiveZ_[31] = 3835.2;
      sensitiveZ_[32] = 3919.1;
      sensitiveZ_[33] = 4003;
      sensitiveZ_[34] = 4086.9;
      sensitiveZ_[35] = 4170.8;
      sensitiveZ_[36] = 4254.7;
      sensitiveZ_[37] = 4338.6;
      sensitiveZ_[38] = 4422.5;
      sensitiveZ_[39] = 4506.4;
      sensitiveZ_[40] = 4590.3;
      etaBoundary_.resize(indices_[6],0);
      for (unsigned iL(0); iL<8; ++iL){
	etaBoundary_[iL] = 1.4;
      }
      etaBoundary_[8] = 1.72042;
      etaBoundary_[9] = 1.81718;
      etaBoundary_[10] = 1.82927;
      etaBoundary_[11] = 1.91612;
      etaBoundary_[12] = 2.02287;
      etaBoundary_[13] = 2.09617;
      etaBoundary_[14] = 2.21281;
      etaBoundary_[15] = 2.28463;
      etaBoundary_[16] = 2.30324;
      etaBoundary_[17] = 2.32153;
      etaBoundary_[18] = 2.33949;
      etaBoundary_[19] = 2.35714;
      etaBoundary_[20] = 2.37449;
      etaBoundary_[21] = 2.39155;
      etaBoundary_[22] = 2.40832;
      etaBoundary_[23] = 2.42483;

      etaBoundary_[25] = 1.72042;
      etaBoundary_[26] = 1.81718;
      etaBoundary_[27] = 1.82927;
      etaBoundary_[28] = 1.91612;
      etaBoundary_[29] = 2.02287;
      etaBoundary_[30] = 2.09617;
      etaBoundary_[31] = 2.21281;
      etaBoundary_[32] = 2.28463;
      etaBoundary_[33] = 2.30324;
      etaBoundary_[34] = 2.32153;
      etaBoundary_[35] = 2.33949;
      etaBoundary_[36] = 2.35714;
      etaBoundary_[37] = 2.37449;
      etaBoundary_[38] = 2.39155;
      etaBoundary_[39] = 2.40832;
      etaBoundary_[40] = 2.42483;
    }
    else if (versionNumber == 62){
      indices_[4] = 0;
      indices_[5] = 16;
      indices_[6] = 16;

      sensitiveZ_.resize(indices_[6],0);
      sensitiveZ_[0] = 3040.5;
      sensitiveZ_[1] = 3091.4;
      sensitiveZ_[2] = 3142.3;
      sensitiveZ_[3] = 3193.2;
      sensitiveZ_[4] = 3277.1;
      sensitiveZ_[5] = 3361;
      sensitiveZ_[6] = 3444.9;
      sensitiveZ_[7] = 3528.8;
      sensitiveZ_[8] = 3612.7;
      sensitiveZ_[9] = 3696.6;
      sensitiveZ_[10] = 3780.5;
      sensitiveZ_[11] = 3864.4;
      sensitiveZ_[12] = 3948.3;
      sensitiveZ_[13] = 4032.2;
      sensitiveZ_[14] = 4116.1;
      sensitiveZ_[15] = 4200;

      etaBoundary_.resize(indices_[6],0);
      etaBoundary_[0] = 1.72042;
      etaBoundary_[1] = 1.81718;
      etaBoundary_[2] = 1.82927;
      etaBoundary_[3] = 1.91612;
      etaBoundary_[4] = 2.02287;
      etaBoundary_[5] = 2.09617;
      etaBoundary_[6] = 2.21281;
      etaBoundary_[7] = 2.28463;
      etaBoundary_[8] = 2.30324;
      etaBoundary_[9] = 2.32153;
      etaBoundary_[10] = 2.33949;
      etaBoundary_[11] = 2.35714;
      etaBoundary_[12] = 2.37449;
      etaBoundary_[13] = 2.39155;
      etaBoundary_[14] = 2.40832;
      etaBoundary_[15] = 2.42483;
    }
    else if (versionNumber == 63){
      indices_[0] = 0;
      indices_[1] = 10;
      indices_[2] = 20;
      indices_[3] = 28;
      indices_[4] = 53;
      indices_[5] = 57;
      indices_[6] = 69;
      sensitiveZ_.resize(indices_[6],0);
      sensitiveZ_[0] = 3198;
      sensitiveZ_[1] = 3207.1;
      sensitiveZ_[2] = 3222.4;
      sensitiveZ_[3] = 3231.5;
      sensitiveZ_[4] = 3246.8;
      sensitiveZ_[5] = 3255.9;
      sensitiveZ_[6] = 3271.2;
      sensitiveZ_[7] = 3280.3;
      sensitiveZ_[8] = 3295.6;
      sensitiveZ_[9] = 3304.7;
      sensitiveZ_[10] = 3320;
      sensitiveZ_[11] = 3329.1;
      sensitiveZ_[12] = 3344.4;
      sensitiveZ_[13] = 3353.5;
      sensitiveZ_[14] = 3368.8;
      sensitiveZ_[15] = 3377.9;
      sensitiveZ_[16] = 3393.2;
      sensitiveZ_[17] = 3402.3;
      sensitiveZ_[18] = 3417.6;
      sensitiveZ_[19] = 3426.7;
      sensitiveZ_[20] = 3442;
      sensitiveZ_[21] = 3451.1;
      sensitiveZ_[22] = 3466.4;
      sensitiveZ_[23] = 3475.5;
      sensitiveZ_[24] = 3490.8;
      sensitiveZ_[25] = 3499.9;
      sensitiveZ_[26] = 3515.2;
      sensitiveZ_[27] = 3524.3;
      sensitiveZ_[28] = 3577.4;
      sensitiveZ_[29] = 3626.4;
      sensitiveZ_[30] = 3675.4;
      sensitiveZ_[31] = 3724.4;
      sensitiveZ_[32] = 3773.4;
      sensitiveZ_[33] = 3822.4;
      sensitiveZ_[34] = 3871.4;
      sensitiveZ_[35] = 3920.4;
      sensitiveZ_[36] = 3969.4;
      sensitiveZ_[37] = 4020.3;
      sensitiveZ_[38] = 4071.2;
      sensitiveZ_[39] = 4122.1;
      sensitiveZ_[40] = 4206;
      sensitiveZ_[41] = 4289.9;
      sensitiveZ_[42] = 4373.8;
      sensitiveZ_[43] = 4457.7;
      sensitiveZ_[44] = 4541.6;
      sensitiveZ_[45] = 4625.5;
      sensitiveZ_[46] = 4709.4;
      sensitiveZ_[47] = 4793.3;
      sensitiveZ_[48] = 4877.2;
      sensitiveZ_[49] = 4961.1;
      sensitiveZ_[50] = 5045;
      sensitiveZ_[51] = 5128.9;
      sensitiveZ_[52] = 0;
      sensitiveZ_[53] = 3971.2;
      sensitiveZ_[54] = 4022.1;
      sensitiveZ_[55] = 4073;
      sensitiveZ_[56] = 4123.9;
      sensitiveZ_[57] = 4207.8;
      sensitiveZ_[58] = 4291.7;
      sensitiveZ_[59] = 4375.6;
      sensitiveZ_[60] = 4459.5;
      sensitiveZ_[61] = 4543.4;
      sensitiveZ_[62] = 4627.3;
      sensitiveZ_[63] = 4711.2;
      sensitiveZ_[64] = 4795.1;
      sensitiveZ_[65] = 4879;
      sensitiveZ_[66] = 4962.9;
      sensitiveZ_[67] = 5046.8;
      sensitiveZ_[68] = 5130.7;

      etaBoundary_.resize(indices_[6],0);
      for (unsigned iL(0); iL<36; ++iL){
	etaBoundary_[iL] = 1.4;
      }
      etaBoundary_[36] = 1.72042;
      etaBoundary_[37] = 1.81718;
      etaBoundary_[38] = 1.82927;
      etaBoundary_[39] = 1.91612;
      etaBoundary_[40] = 2.02287;
      etaBoundary_[41] = 2.09617;
      etaBoundary_[42] = 2.21281;
      etaBoundary_[43] = 2.28463;
      etaBoundary_[44] = 2.30324;
      etaBoundary_[45] = 2.32153;
      etaBoundary_[46] = 2.33949;
      etaBoundary_[47] = 2.35714;
      etaBoundary_[48] = 2.37449;
      etaBoundary_[49] = 2.39155;
      etaBoundary_[50] = 2.40832;
      etaBoundary_[51] = 2.42483;

      etaBoundary_[53] = 1.72042;
      etaBoundary_[54] = 1.81718;
      etaBoundary_[55] = 1.82927;
      etaBoundary_[56] = 1.91612;
      etaBoundary_[57] = 2.02287;
      etaBoundary_[58] = 2.09617;
      etaBoundary_[59] = 2.21281;
      etaBoundary_[60] = 2.28463;
      etaBoundary_[61] = 2.30324;
      etaBoundary_[62] = 2.32153;
      etaBoundary_[63] = 2.33949;
      etaBoundary_[64] = 2.35714;
      etaBoundary_[65] = 2.37449;
      etaBoundary_[66] = 2.39155;
      etaBoundary_[67] = 2.40832;
      etaBoundary_[68] = 2.42483;


      for(size_t i=0; i<sensitiveZ_.size(); i++){
        std::cout << "sensitiveZ_[" << i << "] = " << sensitiveZ_[i] << " etaBoundary=" << etaBoundary_[i] << std::endl;
      }
      for(size_t i=0; i<7; i++){
        std::cout << "indices_[" << i << "] = " << indices_[i] << std::endl;
      }


    }

    else if (versionNumber == 80||versionNumber == 81||versionNumber == 82||versionNumber == 83){
      indices_[0] = 0;
      indices_[1] = 10;
      indices_[2] = 20;
      indices_[3] = 28;
      if (versionNumber == 80){
        std::vector<float> sensZ={3198,3207.1,3222.4,3231.5,3246.8,3255.9,3271.2,3280.3,3295.6,3304.7,3320,3329.1,3344.4,3353.5,3368.8,3377.9,3393.2,3402.3,3417.6,3426.7,3442,3451.1,3466.4,3475.5,3490.8,3499.9,3515.2,3524.3,3577.4,3644.4,3711.4,3778.4,3845.4,3912.4,3979.4,4046.4,4115.3,4184.2,4253.1,4322,4390.9,4459.8,4528.7,4597.6,4666.5,4735.4,4804.3,4873.2,4942.1,5011,5079.9,5148.8,0,4048.2,4117.1,4186,4254.9,4323.8,4392.7,4461.6,4530.5,4599.4,4668.3,4737.2,4806.1,4875,4943.9,5012.8,5081.7,5150.6};
        std::vector<float> etaBound={1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.73444,1.83542,2.02128,2.03732,2.10697,2.22009,2.28845,2.30371,2.31874,2.33355,2.34815,2.36255,2.37674,2.39074,2.40455,2.41817,2.43161,0,1.73444,1.83542,2.02128,2.03732,2.10697,2.22009,2.28845,2.30371,2.31874,2.33355,2.34815,2.36255,2.37674,2.39074,2.40455,2.41817,2.43161};
        indices_[5]=57;
        sensitiveZ_.resize(sensZ.size(),0);
        etaBoundary_.resize(sensZ.size(),0);
        for(size_t i=0; i<sensZ.size(); i++){
          sensitiveZ_[i]=sensZ[i];
          if(sensitiveZ_[i]==0) indices_[4]=i;
          if(sensitiveZ_[i]!=0) etaBoundary_[i]=etaBound[i];
        }          
        indices_[6]=sensZ.size();
      }
      if (versionNumber == 81){
        std::vector<float> sensZ={3198,3207.1,3222.4,3231.5,3246.8,3255.9,3271.2,3280.3,3295.6,3304.7,3320,3329.1,3344.4,3353.5,3368.8,3377.9,3393.2,3402.3,3417.6,3426.7,3442,3451.1,3466.4,3475.5,3490.8,3499.9,3515.2,3524.3,3577.4,3651.4,3725.4,3799.4,3873.4,3947.4,4021.4,4097.3,4210.2,4323.1,4436,4548.9,4661.8,4774.7,4887.6,5000.5,5113.4,5226.3,0,4023.2,4099.1,4212,4324.9,4437.8,4550.7,4663.6,4776.5,4889.4,5002.3,5115.2,5228.1};
        std::vector<float> etaBound={1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.72686,1.82953,2.01705,2.09714,2.22034,2.29846,2.3231,2.34716,2.37067,2.39364,2.41611,2.43809,0,1.72686,1.82953,2.01705,2.09714,2.22034,2.29846,2.3231,2.34716,2.37067,2.39364,2.41611,2.43809};
        indices_[5]=49;
        sensitiveZ_.resize(sensZ.size()+1,0);
        etaBoundary_.resize(sensZ.size()+1,0);
        for(size_t i=0; i<sensZ.size(); i++){
          sensitiveZ_[i]=sensZ[i];
          if(sensitiveZ_[i]==0) indices_[4]=i;
          if(sensitiveZ_[i]!=0) etaBoundary_[i]=etaBound[i];
        }          
        indices_[6]=sensZ.size();
      }
      if (versionNumber == 82){
        std::vector<float> sensZ={3198,3207.1,3222.4,3231.5,3246.8,3255.9,3271.2,3280.3,3295.6,3304.7,3320,3329.1,3344.4,3353.5,3368.8,3377.9,3393.2,3402.3,3417.6,3426.7,3442,3451.1,3466.4,3475.5,3490.8,3499.9,3515.2,3524.3,3577.4,3667.4,3757.4,3847.4,3937.4,4027.4,4119.3,4211.2,4303.1,4395,4486.9,4578.8,4670.7,4762.6,4854.5,4946.4,5038.3,5130.2,0,4029.2,4121.1,4213,4304.9,4396.8,4488.7,4580.6,4672.5,4764.4,4856.3,4948.2,5040.1,5132};
        std::vector<float> etaBound={1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.72448,1.83095,2.02222,2.09737,2.21581,2.28936,2.30965,2.32953,2.34904,2.36817,2.38696,2.4054,2.42352,0,1.72448,1.83095,2.02222,2.09737,2.21581,2.28936,2.30965,2.32953,2.34904,2.36817,2.38696,2.4054,2.42352};
        indices_[5]=49;
        sensitiveZ_.resize(sensZ.size()+1,0);
        etaBoundary_.resize(sensZ.size()+1,0);
        for(size_t i=0; i<sensZ.size(); i++){
          sensitiveZ_[i]=sensZ[i];
          if(sensitiveZ_[i]==0) indices_[4]=i;
          if(sensitiveZ_[i]!=0) etaBoundary_[i]=etaBound[i];
        }          
        indices_[6]=sensZ.size();
      }
      if (versionNumber == 83){
        std::vector<float> sensZ={3198,3207.1,3222.4,3231.5,3246.8,3255.9,3271.2,3280.3,3295.6,3304.7,3320,3329.1,3344.4,3353.5,3368.8,3377.9,3393.2,3402.3,3417.6,3426.7,3442,3451.1,3466.4,3475.5,3490.8,3499.9,3515.2,3524.3,3577.4,3616.7,3656,3695.3,3734.6,3773.9,3813.2,3852.5,3891.8,3931.1,3970.4,4011.6,4052.8,4094,4164.8,4235.6,4306.4,4377.2,4448,4518.8,4589.6,4660.4,4731.2,4802,4872.8,4943.6,5014.4,5085.2,5156,0,3972.2,4013.4,4054.6,4095.8,4166.6,4237.4,4308.2,4379,4449.8,4520.6,4591.4,4662.2,4733,4803.8,4874.6,4945.4,5016.2,5087,5157.8};
        std::vector<float> etaBound={1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.72297,1.81742,1.82722,1.83692,2.01627,2.03283,2.10297,2.21656,2.28539,2.30111,2.31659,2.33184,2.34687,2.36168,2.37627,2.39066,2.40485,2.41884,2.43264,0,1.72297,1.81742,1.82722,1.83692,2.01627,2.03283,2.10297,2.21656,2.28539,2.30111,2.31659,2.33184,2.34687,2.36168,2.37627,2.39066,2.40485,2.41884,2.43264};
        indices_[5]=62;
        sensitiveZ_.resize(sensZ.size()+1,0);
        etaBoundary_.resize(sensZ.size()+1,0);
        for(size_t i=0; i<sensZ.size(); i++){
          sensitiveZ_[i]=sensZ[i];
          if(sensitiveZ_[i]==0) indices_[4]=i;
          if(sensitiveZ_[i]!=0) etaBoundary_[i]=etaBound[i];
        }
        indices_[6]=sensZ.size(); 
      }

      for(size_t i=0; i<sensitiveZ_.size(); i++){
        std::cout << "sensitiveZ_[" << i << "] = " << sensitiveZ_[i] << " etaBoundary=" << etaBoundary_[i] << std::endl;
      }
      std::cout << sensitiveZ_.size() << " " << etaBoundary_.size() << std::endl;
      for(size_t i=0; i<7; i++){
        std::cout << "indices_[" << i << "] = " << indices_[i] << std::endl; 
      }


    }
    else if (versionNumber == 110){
      indices_[0] = 0;
      indices_[1] = 4;
      indices_[2] = indices_[1];
      indices_[3] = indices_[1];
      indices_[4] = indices_[1];
      indices_[5] = indices_[1];
      indices_[6] = indices_[1];
    }
    else {
      indices_[0] = 0;
      indices_[1] = versionNumber==24?11:10;
      indices_[2] = versionNumber==24?21:20;
      indices_[3] = versionNumber==24?31:30;
      indices_[4] = versionNumber==24?55:42;
      indices_[5] = versionNumber==24?65:54;
      indices_[6] = versionNumber==24?65:54;
    }
    
  };

  void buildDetector(const unsigned versionNumber,
		     bool concept=true,
		     bool isCaliceHcal=false,
		     bool bypassR=false);

  const HGCSSSubDetector & subDetectorByLayer(const unsigned aLayer);

  unsigned getSection(const unsigned aLayer) const;
  inline unsigned section(const DetectorEnum adet){
    if (enumMap_.find(adet) != enumMap_.end())
      return enumMap_[adet];
    return nSections();
  };

  void addSubdetector(const HGCSSSubDetector & adet);
  
  void finishInitialisation();

  inline bool isMixedLayer(const unsigned versionNumber,const unsigned aLayer){
    if (versionNumber!=63) return false;
    if (aLayer<36) return false;
    else return true;
  };

  inline unsigned nLayers(const unsigned aSection) const{
    return subdets_[aSection].nLayers();
  };

  inline unsigned nLayers(DetectorEnum adet){
    return subdets_[enumMap_[adet]].nLayers();
  };

  const HGCSSSubDetector & subDetectorByEnum(DetectorEnum adet);
  inline const HGCSSSubDetector & subDetectorBySection(const unsigned aSection) const{
    return subdets_[aSection];
  };

  inline unsigned nLayers() const{
    return nLayers_;
  };

  inline unsigned nSections() const{
    return nSections_;
  };

  inline DetectorEnum detType(const unsigned aSection) const{
    return subdets_[aSection].type;
  };
  
  inline DetectorEnum detTypeLayer(const unsigned aLayer) const{
    return subdets_[getSection(aLayer)].type;
  };
  
  inline std::string detName(const unsigned aSection) const{
    return subdets_[aSection].name;
  };


  void reset();

  void printDetector(std::ostream & aOs) const ;

  inline double sensitiveZ(const unsigned layer){
    if (layer<sensitiveZ_.size()) return sensitiveZ_[layer];
    else {
      std::cout << " ERROR! Trying to access layer " << layer << " outside of range: " << sensitiveZ_.size() << " nLayers " << nLayers_ << std::endl;
      exit(1);
    }
  };

  inline double etaBoundary(const unsigned layer){
    if (layer<etaBoundary_.size()) return etaBoundary_[layer];
    else {
      std::cout << " ERROR! Trying to access layer " << layer << " outside of range: " << etaBoundary_.size() << " nLayers " << nLayers_ << std::endl;
      exit(1);
    }
  };

private:
  HGCSSDetector(){
    bypassRadius_ = false;
  };

  ~HGCSSDetector(){
    reset();
  };
  
  std::vector<HGCSSSubDetector> subdets_;
  std::vector<unsigned> indices_;
  std::vector<unsigned> section_;
  std::map<DetectorEnum,unsigned> enumMap_;

  unsigned nLayers_;
  unsigned nSections_;
  bool bypassRadius_;

  std::vector<double> sensitiveZ_;
  std::vector<double> etaBoundary_;

};

HGCSSDetector & theDetector();




#endif
