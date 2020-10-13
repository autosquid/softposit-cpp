#include "softposit.hpp"
#include <iostream>

namespace bloody
{

    using stat_info = Eigen::Matrix<double, 5, 1>;
    std::tuple<std::vector<unsigned long long>, std::vector<double>>
    maxPosRatio(Eigen::MatrixXd assignMat);
    Eigen::MatrixXd sinkhornImp(Eigen::MatrixXd M);
    int numMatches(Eigen::MatrixXd assignMat);

    std::pair<std::tuple<Pose_type, match_type>, bool> softposit(
        const std::vector<point2di_type> &imagePts,
        const std::vector<point3d_type> &worldPts,
        const Param_type &param,
        const Pose_type &initpose,
        std::pair<const CamInfo_type &, bool> maybe_caminfo)
    {

        using namespace std;

        std::vector<stat_info> stats;

        auto alpha = 9.21 * std::pow(param.noiseStd, 2) + 1;
        auto maxDelta = std::sqrt(alpha) / 2; //  每个世界坐标所允许的最大误差

        auto betaFinal = 0.5;   // 最终迭代目标
        auto betaUpdate = 1.05; // 迭代率
        auto epsilon0 = 0.01;   // 用于初始化任务矩阵

        auto maxCount = 1;
        auto minBetaCount = 20;
        auto nbImagePts = imagePts.size();
        auto nbWorldPts = worldPts.size();

        auto minNbPts = std::min(nbImagePts, nbWorldPts);
        auto maxNbPts = nbImagePts + nbWorldPts - minNbPts;

        auto scale = 1.0 / (maxNbPts + 1);

        std::vector<point2d_type> _centeredImage(imagePts.size());

        CamInfo_type caminfo;
        if (maybe_caminfo.second)
            caminfo = maybe_caminfo.first;
        else
        {
            caminfo.focalLength = 1;
            caminfo.center = point2di_type{0, 0};
        }

        std::cout << "init: " << std::endl
                  << initpose.rot << std::endl
                  << initpose.trans << std::endl;
        Eigen::Matrix<double, Eigen::Dynamic, 2> centeredImage;
        centeredImage.resize(imagePts.size(), Eigen::NoChange);

        for (size_t i = 0; i < imagePts.size(); ++i)
        {
            centeredImage.block<1, 2>(i, 0) << (imagePts[i].x() - (double)caminfo.center.x()) / caminfo.focalLength,
                (imagePts[i].y() - (double)caminfo.center.y()) / caminfo.focalLength;
        }

        std::cout << "centered image :" << centeredImage << std::endl;

        // 齐次坐标
        Eigen::Matrix<double, Eigen::Dynamic, 4> homogeneousWorldPts;
        homogeneousWorldPts.resize(worldPts.size(), Eigen::NoChange);

        for (int i = 0; i < worldPts.size(); ++i)
        {
            homogeneousWorldPts.block<1, 4>(i, 0) << worldPts[i][0], worldPts[i][1], worldPts[i][2], 1;
        }
        std::cout << "begin to make world point homogeneous:" << homogeneousWorldPts << std::endl;

        auto pose = initpose;
        Eigen::VectorXd wk = homogeneousWorldPts * Eigen::Vector4d{pose.rot(2, 0) / pose.trans[2], pose.rot(2, 1) / pose.trans[2], pose.rot(2, 2) / pose.trans[2], 1};
        std::cout << "wk" << wk << std::endl;

        Eigen::Vector4d r1T = {pose.rot(0, 0) / pose.trans(2), pose.rot(0, 1) / pose.trans(2), pose.rot(0, 2) / pose.trans(2), pose.trans(0) / pose.trans(2)};
        Eigen::Vector4d r2T = {pose.rot(1, 0) / pose.trans(2), pose.rot(1, 1) / pose.trans(2), pose.rot(1, 2) / pose.trans(2), pose.trans(1) / pose.trans(2)};

        auto betaCount = 0;
        auto poseConverged = 0;
        auto assignConverged = false;
        auto foundPose = 0;
        auto beta = param.beta0;

        Eigen::MatrixXd assignMat = Eigen::MatrixXd::Ones(nbImagePts + 1, nbWorldPts + 1) * (1 + epsilon0);

        Eigen::VectorXd imageOnes = Eigen::MatrixXd::Ones(nbImagePts, 1);

        int debug_loop = 0;
        while (beta < betaFinal && !assignConverged)
        {
            std::cout << "debug loop: " << debug_loop++ << std::endl;

            Eigen::VectorXd projectedU = homogeneousWorldPts * r1T;
            Eigen::VectorXd projectedV = homogeneousWorldPts * r2T;

            Eigen::MatrixXd replicatedProjectedU = imageOnes * projectedU.transpose();
            Eigen::MatrixXd replicatedProjectedV = imageOnes * projectedV.transpose();

            std::cout << "r1T, r2T used:" << std::endl
                      << r1T << std::endl
                      << r2T << std::endl;
            std::cout << "projected uv:" << std::endl
                      << projectedU << std::endl
                      << projectedV << std::endl;
            std::cout << "reprojected uv" << std::endl;
            std::cout << Eigen::MatrixXd(replicatedProjectedU) << std::endl
                      << Eigen::MatrixXd(replicatedProjectedV) << std::endl;

            std::cout << "SOP" << std::endl;
            auto wkxj = centeredImage.col(0) * wk.transpose();
            auto wkyj = centeredImage.col(1) * wk.transpose();

            std::cout << "wkxj, wkyj" << std::endl;
            std::cout << wkxj << std::endl
                      << wkyj << std::endl;

            Eigen::MatrixXd distMat = caminfo.focalLength * caminfo.focalLength *
                                      ((replicatedProjectedU - wkxj).array().square() +
                                       (replicatedProjectedV - wkyj).array().square());

            std::cout << "dist mat:" << std::endl
                      << distMat << std::endl;

            assignMat.block(0, 0, nbImagePts, nbWorldPts) = scale * (-beta * (distMat - alpha * Eigen::MatrixXd::Ones(distMat.rows(), distMat.cols()))).array().exp();
            assignMat.col(nbWorldPts) = scale * Eigen::VectorXd::Ones(nbImagePts + 1);
            assignMat.row(nbImagePts) = scale * Eigen::VectorXd::Ones(nbWorldPts + 1).transpose();
            std::cout << "assign befor sinkhorn:" << std::endl
                      << assignMat << std::endl;

            assignMat = sinkhornImp(assignMat); // My "improved" Sinkhorn.
            //assignMat = sinkhornSlack (assignMat);    // My "improved" Sinkhorn.
            std::cout << "after sinkhorn Slack:" << std::endl
                      << assignMat << std::endl;

            auto numMatchPts = numMatches(assignMat);
            std::cout << "num matches: " << numMatchPts << std::endl;

            auto sumNonslack = assignMat.block(0, 0, nbImagePts, nbWorldPts).sum();
            std::cout << "sum non slack: " << sumNonslack << std::endl;

            Eigen::VectorXd summedByColAssign = assignMat.block(0, 0, nbImagePts, nbWorldPts).colwise().sum();
            Eigen::Matrix4d sumSkSkT = Eigen::Matrix4d::Zero(4, 4);

            for (auto k = 0; k < nbWorldPts; ++k)
            {
                sumSkSkT = sumSkSkT + summedByColAssign(k) * homogeneousWorldPts.row(k).transpose() * homogeneousWorldPts.row(k);
            }

            std::cout << "check ill-condition" << std::endl;
            Eigen::JacobiSVD<Eigen::Matrix4d> svd_sumSKSkT(sumSkSkT);
            double cond = svd_sumSKSkT.singularValues()(0) / svd_sumSKSkT.singularValues()(svd_sumSKSkT.singularValues().size() - 1);
            if (cond > 1e10)
            {
                std::cout << "sumSkSkT is ill-conditioned, termininating search." << std::endl;
                return std::make_pair(std::make_tuple(pose, match_type()), false);
            }

            Eigen::Matrix4d objectMat = sumSkSkT.inverse(); // Inv(L), a 4x4 matrix.
            poseConverged = 0;                              // Initialize for POSIT loop.
            auto pose_iter_count = 0;

            // Save the previouse pose vectors for convergence checks.
            Eigen::Vector4d r1Tprev = r1T;
            Eigen::Vector4d r2Tprev = r2T;

            double Tx, Ty, Tz;
            Eigen::Vector3d r1, r2, r3;
            double delta;

        std:
            cout << "begin converge loop" << std::endl;
            while (poseConverged == false & pose_iter_count < maxCount)
            {

                Eigen::Vector4d weightedUi = Eigen::Vector4d::Zero(4);
                Eigen::Vector4d weightedVi = Eigen::Vector4d::Zero(4);

                for (int j = 0; j < nbImagePts; ++j)
                {
                    for (int k = 0; k < nbWorldPts; ++k)
                    {
                        weightedUi = weightedUi + assignMat(j, k) * wk(k) * centeredImage(j, 0) * homogeneousWorldPts.row(k).transpose();
                        weightedVi = weightedVi + assignMat(j, k) * wk(k) * centeredImage(j, 1) * homogeneousWorldPts.row(k).transpose();
                    }
                }

                r1T = objectMat * weightedUi;
                r2T = objectMat * weightedVi;

                Eigen::Vector2d s;
                Eigen::Matrix<double, 3, 2> X;

                X.col(0) = r1T.head(3);
                X.col(1) = r2T.head(3);

                std::cout << "svd" << std::endl;
                Eigen::JacobiSVD<Eigen::Matrix<double, 3, 2>> svd(X, Eigen::ComputeFullU | Eigen::ComputeFullV);

                Eigen::Matrix<double, 3, 2> II;
                II << 1, 0, 0, 1, 0, 0;
                s = svd.singularValues();
                Eigen::MatrixXd A = svd.matrixU() * II * svd.matrixV().transpose();

                r1 = A.col(0);
                r2 = A.col(1);
                r3 = r1.cross(r2);

                Tz = 2 / (s(0) + s(1));
                Tx = r1T(3) * Tz;
                Ty = r2T(3) * Tz;
                auto r3T = Eigen::Vector4d{r3[0], r3[1], r3[2], Tz};

                std::cout << "svd" << std::endl
                          << A << r1 << r2 << r3 << Tx << Ty << Tz << std::endl;

                r1T = Eigen::Vector4d{r1[0], r1[1], r1[2], Tx} / Tz;
                r2T = Eigen::Vector4d{r2[0], r2[1], r2[2], Ty} / Tz;
                std::cout << "r1T, r2T update: " << r1T << std::endl
                          << r2T << std::endl;

                wk = homogeneousWorldPts * r3T / Tz;

                std::cout << "delta" << std::endl;
                delta = std::sqrt((assignMat.block(0, 0, nbImagePts, nbWorldPts).cwiseProduct(distMat) / (double)nbWorldPts).sum());
                poseConverged = delta < maxDelta;

                std::cout << "pose converged:" << poseConverged << std::endl;

                std::cout << "generate trace" << std::endl;

                stat_info trace;
                trace << beta, delta, double(numMatchPts) / nbWorldPts,
                    double(sumNonslack) / nbWorldPts,
                    ((r1T - r1Tprev).array().square() + (r2T - r2Tprev).array().square()).sum();

                std::cout << "keep log" << std::endl;
                stats.push_back(trace);

                pose_iter_count = pose_iter_count + 1;
            }

            beta = betaUpdate * beta;
            betaCount = betaCount + 1;
            assignConverged = poseConverged && betaCount > minBetaCount;

            pose.trans = Eigen::Vector3d{Tx, Ty, Tz};
            pose.rot.row(0) = r1.transpose();
            pose.rot.row(1) = r2.transpose();
            pose.rot.row(2) = r3.transpose();

            foundPose = (delta < maxDelta && betaCount > minBetaCount);

            std::cout << "updated pose:" << std::endl
                      << pose.rot << std::endl
                      << pose.trans << std::endl;
            //%% Log:
            std::cout << "converge loop ends" << std::endl;
            std::cout << "pose found: " << foundPose
                      << "delta exit: " << (delta < maxDelta)
                      << ", count exit:" << (betaCount > minBetaCount) << std::endl;
        }

        std::cout << "pose converged:" << std::endl
                  << pose.rot << pose.trans << std::endl;
        return std::make_pair(std::make_tuple(pose, match_type()), true);
    }

    int numMatches(Eigen::MatrixXd assignMat)
    {
        int num = 0;

        auto nimgpnts = assignMat.rows() - 1;
        auto nmodpnts = assignMat.cols() - 1;

        for (int k = 0; k < nmodpnts; ++k)
        {

            unsigned long long imax;
            auto vmax = assignMat.col(k).maxCoeff(&imax);

            if (imax == assignMat.rows() - 1)
                continue; // Slack value is maximum in this column.

            std::vector<unsigned long long> other_cols;
            other_cols.reserve(assignMat.cols() - 1);
            for (int i = 0; i < assignMat.cols(); ++i)
                if (i != k)
                    other_cols.push_back(i);
            unsigned long long kmax;
            assignMat.row(imax).maxCoeff(&kmax);
            if (kmax == k)
            {
                num = num + 1; // This value is maximal in its row & column.
            }
        }
        return num;
    }

    Eigen::MatrixXd sinkhornImp(Eigen::MatrixXd M)
    {
        Eigen::MatrixXd normalizedMat;

        auto iMaxIterSinkhorn = 60; // In PAMI paper
        auto fEpsilon2 = 0.001;     // Used in termination of Sinkhorn Loop.

        auto iNumSinkIter = 0;
        auto nbRows = M.rows();
        auto nbCols = M.cols();

        auto fMdiffSum = fEpsilon2 + 1; // Set "difference" from last M matrix above the loop termination threshold

        // Get the positions and ratios to slack of the nonslack elements that are maximal in their row and column.
        // Eigen::Matrix<unsigned long long, 3, 3>
        std::vector<unsigned long long> posmax;
        std::vector<double> ratios;
        std::tie(posmax, ratios) = maxPosRatio(M);
        // std::cout<<"postmax, rations"<<posmax<<std::endl<<ratios<std::endl;

        while (fabs(fMdiffSum) > fEpsilon2 && iNumSinkIter < iMaxIterSinkhorn)
        {
            auto Mprev = M; // Save M from previous iteration to test for loop termination

            // Col normalization (except outlier row - do not normalize col slacks
            // against each other)
            Eigen::VectorXd McolSums = M.colwise().sum(); // Row vector. 按列求和
            McolSums(nbCols - 1) = 1;                     // Don't normalize slack col terms against each other.

            Eigen::MatrixXd McolSumsRep = Eigen::VectorXd::Ones(nbRows) * McolSums.transpose();
            M = M.cwiseQuotient(McolSumsRep);

            // Fix values in the slack column.
            for (auto i = 0; i < posmax.size(); i += 2)
            {
                M(posmax[i], nbCols - 1) = ratios[i] * M(posmax[i], posmax[i + 1]);
            }
            // Row normalization (except outlier row - do not normalize col slacks against each other)
            Eigen::VectorXd MrowSums = M.rowwise().sum(); // Column vector.
            MrowSums(nbRows - 1) = 1;                     // Don't normalize slack row terms against each other.

            auto MrowSumsRep = MrowSums * Eigen::VectorXd::Ones(nbCols).transpose();
            M = M.cwiseQuotient(MrowSumsRep);

            // Fix values in the slack row.
            for (auto i = 0; i < posmax.size(); i += 2)
            {
                M(nbRows - 1, posmax[i + 1]) = ratios[i + 1] * M(posmax[i], posmax[i + 1]);
            }
            iNumSinkIter = iNumSinkIter + 1;
            fMdiffSum = (M - Mprev).cwiseAbs().sum();
        }
        normalizedMat = M;

        return normalizedMat;
    }

    std::tuple<std::vector<unsigned long long>, std::vector<double>>
    maxPosRatio(Eigen::MatrixXd assignMat)
    {
        //    function [pos, ratios] = maxPosRatio(assignMat)
        std::vector<unsigned long long> pos_;
        std::vector<double> ratios_;

        auto nrows = assignMat.rows(); //size(assignMat,1);
        auto ncols = assignMat.cols(); //size(assignMat,2);
        auto nimgpnts = nrows - 1;
        auto nmodpnts = ncols - 1;
        // 按列遍历
        for (auto k = 0u; k < nmodpnts; ++k)
        {
            unsigned long long imax, __t;
            auto vmax = assignMat.col(k).maxCoeff(&imax);

            if (imax == nrows - 1)
                continue; // Slack value is maximum in this column.

            // 检查列中的最大值在其行中是否为最大值.
            // 在其行中是最大值
            unsigned long long kmax;
            auto _t = assignMat.row(imax).maxCoeff(&kmax);
            if (kmax == k)
            {
                pos_.push_back(imax);
                pos_.push_back(k);
                // Compute the ratios to row and column slack values.
                auto rr = assignMat(imax, ncols - 1) / assignMat(imax, k);
                auto cr = assignMat(nrows - 1, k) / assignMat(imax, k);
                ratios_.push_back(rr);
                ratios_.push_back(cr);
            }
        }
        return std::make_tuple(pos_, ratios_);
    }

    void softposit(float *rot, float *trans, int *foundPose, int *_imagePts, float *_worldPts, int nbImagePts, int nbWorldPts, float beta0, float noiseStd, float *initRot, float *initTrans, float *focalLength, int *center)
    {

        std::vector<point2di_type> imagePts;
        std::vector<point3d_type> worldPts;

        for (uint i = 0u; i < nbImagePts; ++i)
        {
            imagePts.push_back(point2di_type{_imagePts[i * 2], _imagePts[i * 2 + 1]});
        }
        for (uint i = 0u; i < nbWorldPts; ++i)
        {
            worldPts.push_back(point3d_type{_worldPts[i * 3], _worldPts[3 * i + 1], _worldPts[3 * i + 2]});
        }

        Pose_type initpose;
        for (int i = 0, k = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                initpose.rot(i, j) = initRot[k++];
            }
        }
        initpose.trans = point3d_type{initTrans[0], initTrans[1], initTrans[2]};

        Param_type param{beta0, noiseStd};
        if (focalLength)
        {
            CamInfo_type caminfo{*focalLength, point2di_type{*center, *(center + 1)}};
            auto maybe_result = softposit(imagePts, worldPts, param, initpose, std::make_pair(caminfo, true));
        }
        else
        {
            auto maybe_result = softposit(imagePts, worldPts, param, initpose, std::make_pair(CamInfo_type(), false));
        }
    }
} // namespace bloody
