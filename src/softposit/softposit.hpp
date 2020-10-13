#ifndef uu_recreate_wheel_imp_SOFTPOSIT_H
#define uu_recreate_wheel_imp_SOFTPOSIT_H

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

#include <tuple>
#include <vector>
#include <utility>

namespace bloody
{
    // c warp.
    void softposit(float *rot, float *trans, int *foundPose, int *_imagePts, float *_worldPts, int nbImagePts, int nbWorldPts, float beta0, float noiseStd, float *initRot, float *initTrans, float *focalLength, int *center);

    // c++ imp.
    using vec3_type = Eigen::Vector3d;
    using mat3_type = Eigen::Matrix3d;

    using point2di_type = Eigen::Vector2i;
    using point2d_type = Eigen::Vector3d;
    using point3d_type = Eigen::Vector3d;

    using match_type = std::vector<std::pair<int, int>>;

    struct Param_type
    {
        double beta0;
        double noiseStd;
    };

    struct Pose_type
    {
        mat3_type rot;
        vec3_type trans;
    };

    struct CamInfo_type
    {
        double focalLength;
        point2di_type center;
    };

    // at least it will find something
    std::pair<std::tuple<Pose_type, match_type>, bool> softposit(
        const std::vector<point2di_type> &imagePts,
        const std::vector<point3d_type> &worldPts,
        const Param_type &param,
        const Pose_type &initpose,
        std::pair<const CamInfo_type &, bool> maybe_caminfo);

} // namespace bloody

#endif /* uu_recreate_wheel_imp_SOFTPOSIT_H */
