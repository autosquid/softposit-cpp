#include <iostream>
#include <string>
#include <vector>

#include "softposit/softposit.hpp"

int main()
{
  std::string name;
  std::cout<<"hello, "<<name<<std::endl;

  std::vector<bloody::point2di_type_eigen> imagePts;
  std::vector<bloody::point3d_type_eigen> worldPts;

  std::vector<int> raw_imagePts {
    612,
      117,
      486,
      145,
      567,
      206,
      476,
      234,
      441,
      329};

  std::cout<<"load image points"<<std::endl;
  for (uint i=0u; i<raw_imagePts.size(); i+=2){
    imagePts.push_back(bloody::point2di_type_eigen{raw_imagePts[i], raw_imagePts[i+1]});
  }


  std::vector<double> raw_worldPts{ -3.7500,
      0,
      0.5000,
      7.5000,
      0,
      2.7500,
      -3.0000,
      -5.0000,
      -2.0000,
      3.0000,
      5.0000,
      -2.0000,
      0,
      2.2500,
      -0.7500,
      0,
      -2.2500,
      -0.7500
      };

  std::cout<<"load world points"<<std::endl;
  for (uint i=0u; i<raw_worldPts.size(); i+=3){
    worldPts.push_back(bloody::point3d_type_eigen{raw_worldPts[i], raw_worldPts[i+1], raw_worldPts[i+2]});
  }

  bloody::Param_type param{ 2.0E-4, 10.0};
  bloody::CamInfo_type caminfo{982.1f, bloody::point2di_type_eigen{376, 240}};

  bloody::Pose_type initpose;

  initpose.rot << 0, -1, 0, 1, 0, 0, 0, 0, 1;
  initpose.trans = bloody::point3d_type_eigen{0, 0, 30};

  auto maybe_pose = softposit(
    imagePts,
    worldPts,
    param,
    initpose,
    caminfo
    );

  if (maybe_pose){
    auto pose = std::get<0>(*maybe_pose);
    std::cout<<pose.rot<<std::endl;
    std::cout<<pose.trans<<std::endl;
  }else{
    std::cout<<"failed"<<std::endl;
  }

  return 0;

}
