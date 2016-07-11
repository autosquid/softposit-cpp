#ifndef uu_recreate_wheel_imp_SOFTPOSIT_H
#define uu_recreate_wheel_imp_SOFTPOSIT_H

#include <armadillo>
#include <boost/optional.hpp>

#include <tuple>
#include <vector>
#include <utility>


namespace bloody{
// c warp.
  void softposit(float* rot, float* trans, int* foundPose, int* _imagePts, float* _worldPts, int nbImagePts, int nbWorldPts, float beta0, float noiseStd, float* initRot, float* initTrans, float* focalLength, int* center);

// c++ imp.
  using vec3_type = arma::vec3;
  using mat3_type = arma::mat33;

  using point2di_type = arma::ivec2;
  using point2d_type = arma::vec2;
  using point3d_type = arma::vec3;

  using match_type = std::vector<std::pair<int, int> >;

  struct Param_type{
    double beta0;
    double noiseStd;
  };

  struct Pose_type{
    mat3_type rot;
    vec3_type trans;
  };

  struct CamInfo_type{
    float focalLength;
    point2di_type center;
  };


// at least it will find something
  boost::optional< std::tuple<Pose_type, match_type> > softposit(
    const std::vector<point2di_type>& imagePts,
    const std::vector<point3d_type>& worldPts,
    const Param_type& param,
    const Pose_type& initpose,
    boost::optional<const CamInfo_type&> maybe_caminfo);


}

#endif /* uu_recreate_wheel_imp_SOFTPOSIT_H */
