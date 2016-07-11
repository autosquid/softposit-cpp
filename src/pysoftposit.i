%module pysoftposit
%{
#define SWIG_FILE_WITH_INIT
#define SWIG_PYTHON_EXTRA_NATIVE_CONTAINERS

#include "softposit.hpp"

%}

void softposit(float* rot, float* trans, int* foundPose, int* _imagePts, float* _worldPts, int nbImagePts, int nbWorldPts, float beta0, float noiseStd, float* initRot, float* initTrans, float* focalLength, int* center);
