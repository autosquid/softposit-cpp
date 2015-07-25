%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SoftPOSIT Demonstration Code
%
% Date: January 2001 --> April 2003
%
% AUTHORS:
%     Daniel DeMenthon
%     Center for Automation Research
%     University of Maryland
%     College Park, MD 20742-3275
%
%     Philip David
%     Army Research Laboratory        University of Maryland
%     ATTN: AMSRL-CI-CB               Dept. of Computer Science
%     Adelphi, MD 20783               College Park, MD 20742
%
% Copyright 2003, University of Maryland, Army Research Lab, Daniel DeMenthon and Philip David
%
% This program is available under the terms of the GNU General Public License
% as published by the Free Software Foundation.
% See the GNU General Public License for more details: http://www.gnu.org/copyleft/gpl.html
%
% You are free to use this program for non-commercial purpose.
% If you plan to use this code in commercial applications,
% you need additional licensing:
% please contact daniel@cfar.umd.edu or phild@arl.army.mil
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% softPosit    Compute pose and correspondences using the SoftPOSIT algorithm.
%
%    [rot, trans, assignMat, projWorldPts, foundPose, stats] = ...
%        softPosit(imagePts, imageAdj, worldPts, worldAdj, beta0, noiseStd, ...
%                  initRot, initTrans, focalLength, dispLevel, kickout, center)
%
%
% Given a set of 3D world points and a set of 2D image points, determine the
% rotation and translation (the pose) of the world points that best aligns
% the projected world points with the image points.  Correspondences between
% the world and image points are not known.  The set of image points may
% include spurious and missing data (due to clutter and occlusion).
%
% The registration is computed using the SoftPOSIT algorithm, which is
% an iterative algorithm combining the POSIT and SoftAssign algorithms.
% For a technical description of the algorithm, see [D. DeMenthon and
% P. David, "SoftPOSIT: An Algorithm for Registration of 3D Models to Noisy
% Perspective Images Combining SoftAssign and POSIT," Univ.  of Maryland
% Center for Automation Research Technical Report CAR-TR-970, May 2001].
%
% INPUTS:
%     IMAGEPTS is an M x 2 matrix of the x and y coordinates of a set of
% M (M >= 3) image points.
%     IMAGEADJ is an M x M adjacency matrix for the set of image points.
% This is used for the purpose of drawing pictures, not for computing
% the pose of the world points.  If IMAGEADJ(I,J) is nonzero, then image
% point IMAGEPTS(I) will be connected via a line segment to image point
% IMAGEPTS(J) in any images displayed by this code.  If you don't want
% to see any connecting lines, then set IMAGEADJ = ZEROS(M,M).
%     WORLDPTS is an N x 3 matrix of the x, y, and z coordinates of a set of
% N (N >= 4) world points.
%     WORLDADJ is an N x N adjacency matrix for the set of world points.
% This is used for the purpose of drawing pictures, not for computing the
% pose of the world points.  If WORLDADJ(I,J) is nonzero, then the image of
% world point WORLDPTS(I) will be connected via a line segment to the image
% of world point WORLDPTS(J) in any images displayed by this code.  If you
% don't want to see any connecting lines, then set WORLDADJ = ZEROS(N,N).
%     BETA0 defines the initial fuzziness of the correspondence matrix.
% If we know that our initlal pose is close to the actual pose, then
% beta0 should be close to 0.1 so that the correspondence matrix is not
% completely fuzzy;  otherwise, if the pose is completely unknown, beta0
% should be close to 0.0001 so that all possibilities of correspondence
% are possible in the beginnning.
%     NOISESTD is the standard deviation of the normally distributed noise in
% the x and y coordinates of image points.
%     INITROT is a 3 x 3 matrix giving the initial guess for the rotation of
% the world points.  The transformation from world coordinates (Xw) to camera
% coordinates (Xc) is assumed to be Xc = R*Xw + T where R is a rotation matrix
% and T is a translation vector.
%     INITTRANS is a 3-vector giving the initial guess for the translation of
% the world points.
%     FOCALLENGTH is the focal length of the camera.
%     DISPLEVEL is a value in {0,1,2,3,4,5,6} that determines what amount of
% information this routine produces as it runs.  The following table describes
% what output is produced for various values of DISPLEVEL.
%         -----------------------------------------------------------------
%         Function                                         Display Level
%                                                       0  1  2  3  4  5  6
%         -----------------------------------------------------------------
%         Run silently (display nothing).               X  X  X
%         Show short text messages at all iterations.            X  X  X  X
%         Show long text messages at all iterations.                   X  X
%         Show images at all iterations.                            X  X  X
%         Wait for CR press after each iteration.                         X
%         -----------------------------------------------------------------
%     KICKOUT contains information needed to perform an early termination
% test.  That is, the routine uses the information in KICKOUT to
% determine if finding a correct pose is so unlikely that it should
% just terminate so that it can be restarted from different initial
% conditions.  KICKOUT is a structure with two fields: rthldfbeta and
% numMatchable. kickout.numMatchable is an estimate of the number of
% world points that can be matched to image points.  kickout.rthldfbeta
% is an array of length MAXITERATIONS such that, if on iteration K,
% NMP(K)/kickout.numMatchable is less than kickout.rthldfbeta(K), then the
% iteration terminates and returns with FOUNDPOSE == 0. NMP(K) is the number
% of world points matched to image points on iteration K of softPosit.
% If you don't want to use this kickout test, then set kickout.rthldfbeta
% = zeros(1,MAXITERATIONS).  MAXITERATIONS is the maximum number of
% iterations that softPosit() will perform, and can be calculated from
% BETA0, BETAUPDATE, and BETAFINAL.
%     CENTER is a 2-vector giving the x and y coordinate in the image
% of the optical axis.  This argument is optional.  If given, then all
% image points in IMAGEPTS are normalized with respect to CENTER.  If not
% given, then this routine assumes that the image points have already been
% normalizied with respect to this center point.
%
% OUTPUTS:
%     ROT returns the 3 x 3 rotation matrix of the object.  The computed pose
% of the object is defined by ROT and TRANS.
%     TRANS returns the 3-vector translation of the object.
%     ASSIGNMAT returns the assignments of image points to world points.
% ASSIGNMAT is an (M+1) x (N+1) matrix such that image point J corresponds
% to world point K if ASSIGNMAT(J,K) is maximal in row J and column K.
% M and N are as defined above in the description of the inputs.
%     PROJWORLDPTS returns an N x 2 matrix which gives the image (x and y
% coordinates) of all world points for the pose defined by ROT and TRANS.
% If this pose is correct, then most of the rows of PROJWORLDPTS should
% be close to some row of IMAGEPTS.
%     FOUNDPOSE returns nonzero if the softPOSIT iteration converged
% to a pose; otherwise, FOUNDPOSE returns zero.  Note that the fact that
% softPOSIT converged to some pose does not mean that this pose is correct.
% The calling routine should check the assignments returned in ASSIGNMAT
% to decide if this pose provides for a sufficient number of matches.
%     STATS returns an NUMITERS x 5 matrix containing various information
% about the state of softPOSIT on each iteration.  NUMITERS is the number
% of iterations for which softPOSIT ran.  See the code for the exact values
% returned in this matrix.  This information is probably only useful for
% learning the thresholds to use in the early termination test.
%
% EXAMPLE USAGE:
%     This is an example of softPOSIT registering an image of a cube (in
% which one cube vertex is not visible) to a 3D model of the cube.
% In this example, the internal paramter betaUpdate was set to 1.05.
% The pose of the cube that was used to generate the image points is:
%             Rotation = [ -0.0567   -0.9692    0.2397;
%                          -0.7991    0.1879    0.5710;
%                          -0.5985   -0.1592   -0.7852];
%             Translation = [0.5; -0.1; 7];
% The rotation matrix computed by softPOSIT (see below) does not match the
% above rotation matrix because, for a cube, any of a number of different
% rotations register the image to the model equally well.
%        INPUTS:
%            IMAGEPTS = [ 172.3829  -15.4229;
%                         174.9147 -183.8248;
%                         -28.3942 -147.8052;
%                         243.2142  105.4463;
%                         252.6934  -72.3310;
%                          25.7430  -28.9218;
%                          35.9377  149.1948];
%            IMAGEADJ = [ 1     1     0     1     0     0     0;
%                         1     1     1     0     1     0     0;
%                         0     1     1     0     0     1     0;
%                         1     0     0     1     1     0     1;
%                         0     1     0     1     1     1     0;
%                         0     0     1     0     1     1     1;
%                         0     0     0     1     0     1     1];
%            WORLDPTS = [ -0.5000   -0.5000   -0.5000;
%                          0.5000   -0.5000   -0.5000;
%                          0.5000    0.5000   -0.5000;
%                         -0.5000    0.5000   -0.5000;
%                         -0.5000   -0.5000    0.5000;
%                          0.5000   -0.5000    0.5000;
%                          0.5000    0.5000    0.5000;
%                         -0.5000    0.5000    0.5000];
%            WORLDADJ = [ 1     1     0     1     1     0     0     0;
%                         1     1     1     0     0     1     0     0;
%                         0     1     1     1     0     0     1     0;
%                         1     0     1     1     0     0     0     1;
%                         1     0     0     0     1     1     0     1;
%                         0     1     0     0     1     1     1     0;
%                         0     0     1     0     0     1     1     1;
%                         0     0     0     1     1     0     1     1];
%            BETA0 = 2.0e-04
%            NOISESTD = 0
%            INITROT = [ 0.9149    0.1910   -0.3558;
%                       -0.2254    0.9726   -0.0577;
%                        0.3350    0.1330    0.9328];
%            INITTRANS = [0; 0; 50];
%            FOCALLENGTH = 1500;
%            DISPLEVEL = 5;
%            KICKOUT.numMatchable = 7;
%            KICKOUT.rthldfbeta = zeros(1,200);
%
%        EXECUTION:
%            SoftPOSIT runs for 42 iterations.  The message displayed on the
%            final iteration is the following:
%                betaCount = 42, beta = 0.0015523, delta = 0.86674,
%                poseConverged = 1, numMatchPts = 7, sumNonslack = 5.0256
%            On the final iteration, the projected world points exactly overlay
%            the image points.
%
%        OUTPUTS:
%            ROT = [   0.9692    0.0567   -0.2395;
%                     -0.1879    0.7990   -0.5712;
%                      0.1589    0.5987    0.7851];
%            TRANS = [ 0.5002; -0.1001; 7.0011];
%            ASSIGNMAT = [   0     0     0     0     0     0  0.72    0  0.28;
%                            0     0     0     0     0  0.72     0    0  0.28;
%                            0     0     0     0  0.72     0     0    0  0.28;
%                            0     0  0.72     0     0     0     0    0  0.28;
%                            0  0.72     0     0     0     0     0    0  0.28;
%                         0.72     0     0     0     0     0     0    0  0.28;
%                            0     0     0  0.72     0     0     0    0  0.28;
%                         0.28  0.28  0.28  0.28  0.28  0.28  0.28  1.0  0.11];
%            PROJWORLDPTS = [  25.7389  -28.9020;
%                             252.6681  -72.2992;
%                             243.1986  105.4171;
%                              35.9445  149.1457;
%                             -28.3467 -147.8163;
%                             174.9461 -183.8300;
%                             172.4196  -15.4739;
%                             -14.9407   21.2222];
%            FOUNDPOSE = 1;
%            STATS = [0.00020  53.3      0   0.24   0.00108;
%                     0.00021  55.1      0   0.30   0.00112;
%                     0.00022  55.3      0   0.36   0.00105;
%                     0.00023  55.1   0.12   0.44   0.00106;
%                     0.00024  51.2   0.37   0.50   0.00128;
%                     ...];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [rot, trans, assignMat, projWorldPts, foundPose, stats] = ...
%     softPosit(imagePts, imageAdj, worldPts, worldAdj, beta0, noiseStd, ...
%               initRot, initTrans, focalLength, dispLevel, kickout, center)

function [rot, trans] = ...
    softPosit(imagePts, worldPts, beta0, noiseStd, initRot, initTrans, focalLength, center)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check inputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% msg = nargchk(11, 12, nargin);
% if (~isempty(msg))
%         error(strcat('ERROR: ', msg));
% end
%
% clear msg;
%
% if nargin == 11 % image coordinates are already centered, center is not given
%         center = [0, 0];
% end

%Specify center coordinate
%center = [0, 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stats = [];

% Alpha should be set to the maximum allowed squared distance between an
% image point and the matching projected model point.  This should be a
% function of the noise level in the image.  With normally distributed
% x and y noise of mean zero and standard deviation noiseStd, the squared
% distance between a true 2D point and the measured 2D point has a
% chi-squared distribution with 2 degrees of freedom.  Thus, to ensure
% with probability 0.99 that a measured point is allowed to match to a
% true point, we should take alpha = 9.21*noiseStd^2. See Hartley &
% Zisserman, Multiple View Geometry, p. 103 & p. 550.  A value of 1 is
% added to this to account for possible roundoff errors.
alpha = 9.21*noiseStd^2 + 1;

maxDelta = sqrt(alpha)/2;          % Max allowed error per world point.

betaFinal = 0.5;                   % Terminate iteration when beta == betaFinal.
betaUpdate = 1.05;                 % Update rate on beta.
epsilon0 = 0.01;                   % Used to initialize assignement matrix.

maxCount = 2;
minBetaCount = 20;

nbImagePts = size(imagePts, 1);    % Number of image points (M below).
nbWorldPts = size(worldPts, 1);    % Number of world points (N below).

minNbPts = min([nbImagePts;nbWorldPts]);
maxNbPts = max([nbImagePts;nbWorldPts]);
scale = 1/(maxNbPts + 1);

% Convert to normalized image coordinates. With normalized coordinates, (0,0)
% is the point where the optic axis penetrates the image, and the focal length
% is 1.
centeredImage(:,1) = (imagePts(:,1) - center(1))/focalLength;
centeredImage(:,2) = (imagePts(:,2) - center(2))/focalLength;

ipose.mageOnes = ones(nbImagePts, 1);   % Column M-vector.

% The homogeneous world points -- append a 1 to each 3-vector. An Nx4 matrix.
worldOnes = ones(nbWorldPts, 1);   % Column N-vector.
homogeneousWorldPts = [worldPts, worldOnes];

% Initial rotation and translation as passed into this function.
rot = initRot;
trans = initTrans;
% if ismember(dispLevel,[5,6])
%     disp(' '); disp('Initial transformation:'); rot, trans
% end

% Initialize the depths of all world points based on initial pose.
wk = homogeneousWorldPts * [rot(3,:)/trans(3), 1]';

% Draw a picture of the model on the original image plane.
%projWorldPts = proj3dto2d(worldPts, rot, trans, focalLength, 1, center);
% if ismember(dispLevel,[4,5,6])
%     % plotImages(imagePts, imageAdj, projWorldPts, worldAdj);
%     plot2dPts(imagePts, 'r.-', 'FigNum', 1, 'AdjMat', imageAdj, ...
%               'Name', 'Original image', 'Axis', 300);
%     plot2dPts(projWorldPts, 'b.-', 'Overlay', 'AdjMat', worldAdj, ...
%               'Name', 'Projected model', 'Axis', 'Freeze');
% end

% First two rows of the camera matrices (for both perspective and SOP).  Note:
% the scale factor is s = f/Tz = 1/Tz since f = 1.  These are column 4-vectors.
r1T = [rot(1,:)/trans(3), trans(1)/trans(3)]';
r2T = [rot(2,:)/trans(3), trans(2)/trans(3)]';

betaCount = 0;
poseConverged = 0;
assignConverged = 0;
foundPose = 0;
beta = beta0;

% The assignment matrix includes a slack row and slack column.
assignMat = ones(nbImagePts+1,nbWorldPts+1) + epsilon0;
%keyboard;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deterministic annealing loop.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while beta < betaFinal & ~assignConverged

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Given the current pose and w[i], compute the distance matrix, d[j,k].
    % d[j,k] is the squared distance between the j-th corrected SOP image
    % point and the SOP projection of the k-th scene point.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % SOP (weak perspective projection) of world points.
    projectedU = homogeneousWorldPts * r1T;          % Column N-vector.
    projectedV = homogeneousWorldPts * r2T;          % Column N-vector.

    % MxN matrices with identical rows equal to the X (replicatedProjectedU)
    % and Y (replicatedProjectedV) SOP projections of the world points.  The
    % number of rows M is equal to the number of image points.
    replicatedProjectedU = imageOnes * projectedU';
    replicatedProjectedV = imageOnes * projectedV';

    % For all possible pairings of image to world points, map the
    % perspective image point into the corrected SOP image point using
    % the world point's current estimate of w[i].  Here, j is the index of
    % an image point and k is the index of a world point.  wkxj is an
    % MxN  matrix where the j-th,k-th entry is w[k]*x[j]/f.  wkyj is an
    % MxN matrix where the j-th,k-th entry is w[k]*y[j]/f.
    wkxj = centeredImage(:,1) * wk';
    wkyj = centeredImage(:,2) * wk';

    % distMat[j,k] = d[j,k]^2.  The focalLength^2 term here maps the distances
    % back to the original (unnormalized) image, so distances are in terms
    % of original pixels.
    distMat = focalLength*focalLength*((replicatedProjectedU - wkxj).^2 + ...
                                       (replicatedProjectedV - wkyj).^2);
%     if ismember(dispLevel,[5,6])
%         disp(' '); disp('Distance matrix ='); disp(distMat);
%     end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Use the softAssign algorithm to compute the best assignment of image to
    % world points based on the computed distances.  The use of alpha
    % determines when to favor assignments instead of use of slack.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    assignMat(1:nbImagePts, 1:nbWorldPts) = scale*exp(-beta*(distMat - alpha));
    assignMat(1:nbImagePts+1,nbWorldPts+1) = scale;
    assignMat(nbImagePts+1,1:nbWorldPts+1) = scale;
    % assignMat = sinkhornSlack(assignMat); % Don't normalize slack row and col.

    %assignMat
    %keyboard;

    assignMat = sinkhornImp(assignMat);    % My "improved" Sinkhorn.

    %assignMat
    %keyboard;

%     if ismember(dispLevel,[5,6])
%         disp(' '); disp('Processed assignment matrix ='); disp(assignMat);
%     end

    % About how many matching model points do we currently have?
    numMatchPts = numMatches(assignMat);
    sumNonslack = sum(sum(assignMat(1:nbImagePts,1:nbWorldPts)));

    %keyboard;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Use POSIT to calculate the pose that minimizes the objective function
    % (equation (10), the weighted sum of the differences between corrected
    % image points and SOP's of the corresponding world points), given the
    % current assignment.  The pose parameters and the w[i] are iteratively
    % updated until some convergence criterion is satisfied.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Compute Sum{k=1,N}(m'[k]S[k]S[k]^T) -- this is used in equations (11) and
    % (12) in the paper, and below in the POSIT loop, to compute an optimal
    % pose.  In the paper, L == sumSkSkT.
    %     First, compute the sum of the elements w/in each column,
    % ignoring the slack row and column. This is a row N-vector.
    summedByColAssign = sum(assignMat(1:nbImagePts, 1:nbWorldPts), 1);
    %     Now, compute the 4x4 matrix L.
    sumSkSkT = zeros(4, 4);
    for k = 1:nbWorldPts
        sumSkSkT = sumSkSkT + summedByColAssign(k) * ...
                        homogeneousWorldPts(k,:)' * homogeneousWorldPts(k,:);
    end

    if cond(sumSkSkT) > 1e10
        disp('sumSkSkT is ill-conditioned, termininating search.');
        return;
    end

    objectMat = inv(sumSkSkT);                           % Inv(L), a 4x4 matrix.


    poseConverged = 0;                              % Initialize for POSIT loop.
    count = 0;

    % Save the previouse pose vectors for convergence checks.
    r1Tprev = r1T;
    r2Tprev = r2T;



    % POSIT loop.  We put a cap on number of steps here, so at first it does
    % not converge.  We currently do just one iteration of POSIT.  When the
    % annealing temperature is slow enough, one iteration of POSIT is
    % sufficient since the assigments and pose will converge simultaneously.
     while ~poseConverged & count < maxCount

        % Compute the term in the equation for the optimal pose vectors (M and
        % N in equations (11) and (12) in the paper) that depends on the current
        % assigment and the w[i]. These are Sum{all j,k}(m[j,k]w[k]x[j]S[k]) and
        % Sum{all j,k}(m[j,k]w[k]y[j]S[k]).  Here, (w[k]x[j],w[k]y[j]) is our
        % best estimate of the SOP of the k-th scene point whose homogeneous
        % perspective image coordinate is (x[j],y[j]).
        weightedUi = zeros(4,1);                    % column vector
        weightedVi = zeros(4,1);                    % column vector
        for j = 1:nbImagePts
            for k = 1:nbWorldPts
                weightedUi = weightedUi + assignMat(j,k) * wk(k) * ...
                                 centeredImage(j,1) * homogeneousWorldPts(k,:)';
                weightedVi = weightedVi + assignMat(j,k) * wk(k) * ...
                                 centeredImage(j,2) * homogeneousWorldPts(k,:)';
            end % for k
        end % for j



        % Compute the pose vectors. M = s(R1,Tx) and N = s(R2,Ty) where the
        % scale factor is s = f/Tz, and f = 1.  These are column 4-vectors.
        r1T= objectMat * weightedUi;               % M
        r2T = objectMat * weightedVi;              % N

        % Compute the rotation matrix and translation vector corresponding to
        % the computed pose vectors.
        if 1  % Chang & Tsai calculation of R and T.
            [U, S, V] = svd([r1T(1:3)'; r2T(1:3)']');
            A = U * [1 0; 0 1; 0 0] * V';
            r1 = A(:,1);
            r2 = A(:,2);
            r3 = cross(r1,r2);
            Tz = 2 / (S(1,1) + S(2,2));
            Tx = r1T(4) * Tz;
            Ty = r2T(4) * Tz;
            r3T= [r3; Tz];


        else
            % Standard calculation of R and T.  The rotation matrix may not be
            % orthonormal.  The object must distort when the rotation matrix
            % is not orthonormal.
            r1TSquared = r1T(1)*r1T(1) + r1T(2)*r1T(2)+ r1T(3)*r1T(3);
            r2TSquared = r2T(1)*r2T(1) + r2T(2)*r2T(2)+ r2T(3)*r2T(3);
            Tz = sqrt(2/(r1TSquared+r2TSquared));   % Chang & Tsai's suggestion.
            r1N = r1T*Tz;                   % Column 4-vectors: (R1,Tx).
            r2N = r2T*Tz;                   %                   (R2,Ty).
            r1 = r1N(1:3);                  % Three rows of the rotation matrix.
            r2 = r2N(1:3);
            r3 = cross(r1,r2);
            r3T= [r3; Tz];                  % Column 4-vector: (R3,Tz).
            Tx = r1N(4);
            Ty = r2N(4);
        end


        % Make the pose vectors consistent with the new rotation matrix.
        % These are column 4-vectors.
        r1T = [r1; Tx]/Tz;
        r2T = [r2; Ty]/Tz;

        % Compute updated estimates for the w[i] (the third homogeneous
        % coordinate of the perspective image points).  The improved w[i] are
        % used above to map (correct) the perspective image coordinates into
        % the corresponding scaled orthographic projection image coordinates.
        wk = homogeneousWorldPts * r3T /Tz;

        % Determine if the pose as computed by POSIT has converged.
        % The pose is considered to have converged when the sum of the
        % weighted distances between the corrected SOP image points and
        % the SOP's of the world points is less than some threshold.
        % This distance is in terms of pixels in the original image
        % coordinate frame.
        delta = sqrt(sum(sum(assignMat(1:nbImagePts,1:nbWorldPts) ...
                                        .* distMat))/nbWorldPts);
        poseConverged = delta < maxDelta;

%         if ismember(dispLevel,[3,4,5,6])
%             % Display some information about the state of the iteration.
%             disp(['betaCount = ' num2str(betaCount) ...
%                   ', beta = ' num2str(beta) ...
%                   ', delta = ' num2str(delta) ...
%                   ', poseConverged = ' num2str(poseConverged) ...
%                   ', numMatchPts = ' num2str(numMatchPts) ...
%                   ', sumNonslack = ' num2str(sumNonslack) ]);
%         end

%         if ismember(dispLevel,[4,5,6])
%             pause;
%         end

        % Save some information for use by the calling routine.
        stats(betaCount+1,1) = beta;
        stats(betaCount+1,2) = delta;
        stats(betaCount+1,3) = numMatchPts/nbWorldPts;
        stats(betaCount+1,4) = sumNonslack/nbWorldPts;
        stats(betaCount+1,5) = sum([(r1T-r1Tprev)' (r2T-r2Tprev)'].^2);

        count = count + 1;                         % Number of POSIT iterations.

    end % of POSIT loop

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updates for deterministic annealing and checks for convergence.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Update the "annealing temperature" and determine if the assignments
    % have converged.
    beta = betaUpdate * beta;
    betaCount = betaCount + 1;   % Number of deterministic annealing iterations.
    assignConverged = poseConverged & betaCount > minBetaCount;

    % Form the best estimates for the translation vector and rotation matrix.
    trans = [Tx; Ty; Tz];
    rot = [r1'; r2'; r3'];

%     if ismember(dispLevel,[5,6])
%         disp(' '); disp('Current transformation:'); rot, trans
%     end

    % Has the pose converged?
    foundPose = (delta < maxDelta & betaCount > minBetaCount);


    % Project the model onto the original (unnormalized) image plane.
    %projWorldPts = proj3dto2d(worldPts, rot, trans, focalLength, 1, center);
%     if ismember(dispLevel,[4,5,6])
%         % plotImages(imagePts, imageAdj, projWorldPts, worldAdj);
%         plot2dPts(imagePts, 'r.-', 'FigNum', 1, 'AdjMat', imageAdj, ...
%                   'Name', 'Original image');
%         plot2dPts(projWorldPts, 'b.-', 'Overlay', 'AdjMat', worldAdj, ...
%                   'Name', 'Projected model');
%     end

    % Perform the "early restart" test.
%     r = numMatchPts/kickoutNumMatchable;
%     if r < kickoutRthldfbeta(betaCount)
%         if ismember(dispLevel,[3,4,5,6])
%             disp(['Early restart: r = ' num2str(r) ' < r*(' ...
%                    num2str(betaCount) ') = ' ...
%                    num2str(kickoutRthldfbeta(betaCount))]);
%         end
%         return;
%     end

end % of Deterministic annealing loop

% if poseConverged
%     % Draw a picture showing the initial and final pose of the model.
%     initProjWorldPts = proj3dto2d(worldPts, initRot, initTrans, ...
%                                   focalLength, 1, center);
%     finalProjWorldPts = proj3dto2d(worldPts,rot,trans,focalLength,1,center);
%     plot2dPts(initProjWorldPts, 'g.-', 'FigNum', 5, 'AdjMat', worldAdj, ...
%               'Axis', 'Auto');
%     plot2dPts(finalProjWorldPts, 'b.-', 'Overlay', 'AdjMat', worldAdj);
%     plot2dPts(imagePts, 'r.', 'Overlay', 'AdjMat', imageAdj);
%     axis(axis+[-50 50 -50 50]);                    % Expand the axis a little.
%     set(gca,'XTick',[])
%     set(gca,'YTick',[])
%     print('-depsc','-f5',['poses' num2str(fix(rand*1000)) '.ps']);
%     disp('A picture of the initial and final poses has been saved.');
% end

return






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROJ3DTO2D   Project a set of 3D points onto a 2D image.
%      PTS2D = proj3dto2d(PTS3D, ROT, TRANS, FLENGTH, OBJDIM, [CENTER])
%      returns a Nx2 (or 2xN, see below) matrix of the projections of the
%      3D object points PTS3D onto the image plane at focal length FLENGTH.
%      The 3D object points are mapped into the camera coordinate system
%      according to the transformation "Xcamera = ROT*Xobject + TRANS".
%      OBJDIM gives the dimension (1 or 2) in PTS3D along which object
%      points are indexed. Thus, the image matrix PTS2D will be Nx2
%      if OBJDIM == 1 and 2xN if OBJDIM == 2.  The optional argument
%      CENTER is a 2-vector which is principle point of the image plane;
%      if not specified, the center is assumed to be (0,0).
%
%    AUTHOR:
%        Philip David
%        Army Research Laboratory        University of Maryland
%        ATTN: AMSRL-CI-CB               Dept. of Computer Science
%        Adelphi, MD 20783               College Park, MD 20742
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function pts2d = proj3dto2d(pts3d, rot, trans, flength, objdim, varargin)
%
% % Check the inputs.
% if objdim ~= 1 & objdim ~= 2
%     error(['OBJDIM must be 1 or 2.  You gave ' num2str(OBJDIM)]);
% elseif size(pts3d,bitcmp(objdim,2)) ~= 3
%     error(['PTS3D does not have size 3 along dimension ' ...
%           num2str(bitcmp(objdim,2))]);
% end
%
% % Get the principle point of the image.
% if length(varargin) >= 1
%     center = varargin{1};
%     if size(center,2) ~= 1
%         center = center';              % Make sure center is a column vector.
%     end
% else
%     center = [0;0];
% end
%
% % Make sure the 3D point matrix is 3xN.
% if objdim == 1
%     pts3d = pts3d';
% end
%
% numpts = size(pts3d,2);       % Number of 3D points.
%
% % Map the world points into camera points.
% campts = rot*pts3d + trans(:,ones(1,numpts));          % 3D camera points.
% campts(3,find(campts(3,:) < 1e-20)) = 1e-20;           % Don't divide by zero!
%
% % Project the camera points into the image.
% pts2d = flength * campts(1:2,:) ./ campts([3;3],:);    % 2D image points.
% pts2d = pts2d + center(:,ones(1,numpts));              % Shift principle point.
%
% % Put the 2D point matrix into the same form as the 3D point matrix.
% if objdim == 1
%     pts2d = pts2d';
% end
%
% return





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DISP2DPTS   Plot a set of 2D points in a figure.
%
%    plot2dPts(PNTPOS, LINESPEC, [OPTIONAL_ARGS]) plots the 2D points
%    in PNTPOS using the line specification LINESPEC.  PNTPOS is an Nx2
%    or a 2xN matrix of the X and Y coordinates of the points to plot.
%    (If N == 2, then each row of PNTPOS is assumed to be a different
%    2D point.)  LINESPEC is a character string of the type passed to
%    PLOT() that specifies the line color, marker type, and line type.
%    Most optional arguments are remembered from one call of this routine
%    to the next (and are associated with particular figure numbers), and
%    therefore don't need to be repeated from one call to the next.
%
%    Various optional arguments may also appear in any order:
%        'Title', TITLE      -- String giving title of the figure.
%        'Name', NAME        -- String giving name of the current data set.
%        'AdjMat', ADJMAT    -- NxN adjacency matrix for the data.
%        'FigNum', FIGNUM    -- Figure number. FIGNUM may be a single integer
%                               or may be [FN, NR, NC, POS]. In the later case,
%                               FN is the figure number, the figure has
%                               NRxNC axes, and POS is the position for
%                               the current plot.
%        'New'               -- Start a new figure?
%        'Overlay'           -- Overlay current points ontop of previous points?
%        'LineWidth', LWIDTH -- Line width.
%        'LabelPts'          -- Label the points?
%        'Axis', AXIS        -- AXIS may be a single number, the vector
%                               [XMIN XMAX YMIN YMAX], the string 'AUTO',
%                               or the string 'FREEZE'.  For the case
%                               of a single number, EACH axis displays
%                               the range -AXIS to +AXIS.  For the
%                               case of a vector, the axis are set to
%                               [XMIN XMAX YMIN YMAX].  When AXIS is
%                               'AUTO', the axis are automatically set.
%                               When AXIS is 'FREEZE', the axis are frozen
%                               after plotting the current data.
%
%    AUTHOR:
%        Philip David
%        Army Research Laboratory        University of Maryland
%        ATTN: AMSRL-CI-CB               Dept. of Computer Science
%        Adelphi, MD 20783               College Park, MD 20742
%
%    HISTORY:
%        14 July 2001 - P. David - Created.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function plot2dPts(pntpos, linespec, varargin)
%
% persistent firstcall;
% persistent fignum;
% persistent figpos;
% persistent numpntsets;
% persistent pntnames;
% persistent numlabeledsets;
% persistent figtitle;
% persistent fignrows;
% persistent figncols;
% persistent linewidth;
% persistent axismode;           % 0 = No change, 1 = Auto, 2 = Manual
% persistent axisxmin;
% persistent axisxmax;
% persistent axisymin;
% persistent axisymax;
% persistent showlegend;
%
% % Initialize the default parameters.
% if isempty(firstcall)
%     firstcall = 0;
%     fignum = 1;
%     figpos = 1;
%     for k = 1:100
%         pntnames{k} = [];
%         for j = 1:20
%             numlabeledsets{k,j} = 0;
%             numpntsets{k,j} = 0;
%             showlegend{k,j} = 0;
%         end
%         figtitle{k} = [];
%         fignrows{k} = 1;
%         figncols{k} = 1;
%         linewidth{k} = 1;
%         axismode{k} = 1;
%         axisxmin{k} = -500;
%         axisxmax{k} = 500;
%         axisymin{k} = -500;
%         axisymax{k} = 500;
%     end
% end
%
% % Make sure the data matrix PNTPOS is an Nx2 matrix.
% [nrows, ncols] = size(pntpos);
% if nrows == 0 & ncols == 0
%     % We allow an empty set of points.
% elseif nrows ~= 2 & ncols ~= 2
%     error('PNTPOS must be Nx2 or 2xN');
% elseif nrows == 2 & ncols == 2
%     disp('plot2dPts: 2x2 point matrix, assuming each row is a 2D point.');
%     pause;
% elseif nrows == 2
%     % Make the 2D point matrix is Nx2;
%     pntpos = pntpos';
%     tmp = ncols;
%     ncols = nrows;
%     nrows = ncols;
% end
%
% % Set default arguments.
% curname = '';                  % Data set is not named.
% pntadj = eye(nrows);           % No point is connected to any other.
% newfig = 0;                    % Use previous figure number.
% overlayfig = 0;                % Erase previous plot in this figure.
% labelpts = 0;                  % Do not label the current set of points.
%
% % Get optional arguments.
% argnum = 1;
% while argnum <= length(varargin)
%     switch lower(varargin{argnum})
%     case 'title'
%         argnum = argnum + 1;
%         figtitle{fignum} = varargin{argnum};
%     case 'name'
%         argnum = argnum + 1;
%         curname = varargin{argnum};
%     case 'adjmat'
%         argnum = argnum + 1;
%         pntadj = varargin{argnum};
%     case 'fignum'
%         argnum = argnum + 1;
%         fignumdata = varargin{argnum};
%         if length(fignumdata) == 1
%             % Make a figure with one axis.
%             fignum = fignumdata;
%             fignrows{fignum} = 1;
%             figncols{fignum} = 1;
%             figpos = 1;
%         elseif length(fignumdata) == 4
%             % Make a figure with an MxN matrix of axes.
%             fignum = fignumdata(1);
%             fignrows{fignum} = fignumdata(2);
%             figncols{fignum} = fignumdata(3);
%             figpos = fignumdata(4);
%         else
%             error(['Invalid FIGNUM argument: ' mat2str(fignum)]);
%         end
%     case 'new'
%         newfig = 1;
%     case 'overlay'
%         overlayfig = 1;
%     case 'axis'
%         argnum = argnum + 1;
%         axisdata = varargin{argnum};
%         if strcmpi(axisdata, 'AUTO')
%             % Let the plot function pick a good axis.
%             axismode{fignum} = 1;
%         elseif strcmpi(axisdata, 'FREEZE')
%             % Freeze the axis after plotting the current data.
%             axismode{fignum} = 0;
%         elseif length(axisdata) == 1
%             % Display the same zero-centered range on both axes.
%             axismode{fignum} = 2;
%             axisxmin{fignum} = -axisdata;
%             axisxmax{fignum} = axisdata;
%             axisymin{fignum} = -axisdata;
%             axisymax{fignum} = axisdata;
%         elseif length(axisdata) == 4
%             % Get the axis info from the vector [XMIN XMAX YMIN YMAX].
%             axismode{fignum} = 2;
%             axisxmin{fignum} = axisdata(1);
%             axisxmax{fignum} = axisdata(2);
%             axisymin{fignum} = axisdata(3);
%             axisymax{fignum} = axisdata(4);
%         else
%             error(['Invalid AXIS argument: ' mat2str(axisdata)]);
%         end
%     case 'linewidth'
%         argnum = argnum + 1;
%         linewidth{fignum} = varargin{argnum};
%     case 'labelpts'
%         labelpts = 1;
%     otherwise
%         error(['Unknown argument label: ''' varargin{argnum} '''']);
%     end
%     argnum = argnum + 1;
% end
%
% figure(fignum);
%
% % Start a new figure?
% if newfig
%     pntnames{fignum} = [];
%     for j = 1:20
%         numlabeledsets{fignum,j} = 0;
%         numpntsets{fignum,j} = 0;
%     end
%     clf;
% end
%
% if overlayfig & numpntsets{fignum,figpos} > 0
%     % If there was a previous plot for this axis, then overlay the
%     % current plot on top of it.
%     numpntsets{fignum,figpos} = numpntsets{fignum,figpos} + 1;
%     hold on;
% else
%     % Draw the first plot for this axis.
%     pntnames{fignum,figpos} = [];
%     numpntsets{fignum,figpos} = 1;
%     numlabeledsets{fignum,figpos} = 0;
%     showlegend{fignum,figpos} = 0;
%     hold off;
% end
%
% subplot(fignrows{fignum}, figncols{fignum}, figpos);
% set(fignum,'DoubleBuffer','on');
%
% % Plot the data if there is any data to plot.
% if nrows > 0, gplot(pntadj, pntpos, linespec); end
%
% % Set the range of data to be displayed.
% axis equal;
% if axismode{fignum} == 1
%     axis auto;
% elseif axismode{fignum} == 2
%     axis([axisxmin{fignum} axisxmax{fignum} axisymin{fignum} axisymax{fignum}]);
% end
%
% % Update the legend.
% %     In the past, putting a legend on the image causes some of the
% %     markers to disappear!  See below for the makeshift legend.
% pntnames{fignum,figpos,numpntsets{fignum,figpos}} = curname;
% if ~isempty(curname) | showlegend{fignum,figpos}
%     showlegend{fignum,figpos} = 1;
%     legend(pntnames{fignum,figpos,1:numpntsets{fignum,figpos}},0);
% end
%
% % Tag the current data according to it's FIGNUM, FIGPOS, and PNTSETNUM
% h = findobj('Type','Line','Color',linespec(1),'Tag','');
% set(h,'Tag',mat2str([fignum figpos numpntsets{fignum,figpos}]));
%
% % Set the size of all markers and lines that are the current color.
% markersize = 8*ceil(linewidth{fignum});
% set(h, 'MarkerSize', markersize);
% set(h, 'LineWidth', ceil(linewidth{fignum}));
%
% % Get the width and height of the plot.
% axisrange = axis;                           % [xmin xmax ymin ymax]
% width = axisrange(2) - axisrange(1);
% height = axisrange(4) - axisrange(3);
% set(gca,'units','points');
% axsize = get(gca,'position');           % [left bottom width height] in points.
% markratio = markersize./axsize(3:4);    % marker-size-to-axis-size ratio.
%
% % Label the points.
% if labelpts
%     switch mod(numlabeledsets{fignum,figpos},4)
%     case 0
%         % halign = 'left'; valign = 'top';
%         % dx = 1; dy = -1;
%         halign = 'left'; valign = 'middle';
%         dx = 1; dy = 0;
%     case 1
%         % halign = 'right'; valign = 'bottom';
%         % dx = -1; dy = 1;
%         halign = 'right'; valign = 'middle';
%         dx = -1; dy = 0;
%     case 2
%         % halign = 'right'; valign = 'top';
%         % dx = -1; dy = -1;
%         halign = 'center'; valign = 'top';
%         dx = 0; dy = -1;
%     case 3
%         % halign = 'left'; valign = 'bottom';
%         % dx = 1; dy = -1;
%         halign = 'center'; valign = 'bottom';
%         dx = 0; dy = 1;
%     end
%     dx = 0.5*dx*markratio(1)*width;
%     dy = 0; % 0.4*dy*markratio(2)*height;
%     for k = 1:nrows
%         h = text(pntpos(k,1)+dx, pntpos(k,2)+dy, int2str(k), ...
%                  'VerticalAlignment', valign, 'HorizontalAlignment', halign, ...
%                  'FontSize', 10, 'Color', linespec(1));
%     end
%     numlabeledsets{fignum,figpos} = numlabeledsets{fignum,figpos} + 1;
% end
%
% % Create a makeshift legend until the regular legend works.
% % if ~isempty(curname)
% %     text(axisrange(1)+0.05*width, ...
% %          axisrange(4)-(numpntsets{fignum,figpos})*0.05*height, ...
% %          curname, 'FontSize', 10, 'Color', linespec(1));
% % end
%
% % Display a figure title.
% if ~isempty(figtitle{fignum})
%     title(figtitle{fignum});
% end
%
% drawnow;
%
% if axismode{fignum} == 0
%     % Freeze the axis at their current settings.
%     axisxmin{fignum} = axisrange(1);
%     axisxmax{fignum} = axisrange(2);
%     axisymin{fignum} = axisrange(3);
%     axisymax{fignum} = axisrange(4);
%     axismode{fignum} = 2;
% end
%
% return;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Normalize across rows and columns to find an assignment matrix (a
% doubly stochastic matrix).
%
% This version treats the slack row and column differently from other rows
% and columns: the slack values are not normalized with respect to other
% slack values, only with respect to the nonslack values.  This may work
% better than the original Sinkhorn algorithm which treats all rows and
% columns identically.  This is true primarily when there needs to be more
% than one assignment to a slack row or column.  I.e., when there are two
% or more missing image points or model points.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function normalizedMat = sinkhornSlack(M)

iMaxIterSinkhorn=60;  % In PAMI paper
fEpsilon2 = 0.001; % Used in termination of Sinkhorn Loop.

iNumSinkIter = 0;
[nbRows, nbCols] = size(M);

fMdiffSum = fEpsilon2 + 1;  % Set "difference" from last M matrix above
                            % the loop termination threshold.

while(abs(fMdiffSum) > fEpsilon2 & iNumSinkIter < iMaxIterSinkhorn)
    Mprev = M;  % Save M from previous iteration to test for loop termination

    % Col normalization (except outlier row - do not normalize col slacks
    % against each other)
    McolSums = sum(M, 1); % Row vector.
    McolSums(nbCols) = 1; % Don't normalize slack col terms against each other.
    McolSumsRep = ones(nbRows,1) * McolSums ;
    M = M ./ McolSumsRep;

    % Row normalization (except outlier row - do not normalize col slacks
    % against each other)
    MrowSums = sum(M, 2); % Column vector.
    MrowSums(nbRows) = 1; % Don't normalize slack row terms against each other.
    MrowSumsRep = MrowSums * ones(1, nbCols);
    M = M ./ MrowSumsRep;

    iNumSinkIter=iNumSinkIter+1;
    fMdiffSum=sum(abs(M(:)-Mprev(:)));
end

normalizedMat = M;

return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% normalizedMat = sinkhornImp(M)
%
% Apply an improved Sinkhorn algorithm to map matrix M to a doubly
% stochastic matrix.
%
% The Sinkhorn algorithm modified for slack rows and columns treats the
% slack row and column differently from other rows and columns: the slack
% values are not normalized with respect to other slack values, only with
% respect to the nonslack values.  This may work better than the original
% Sinkhorn algorithm which treats all rows and columns identically.
% This is true primarily when there needs to be more than one assignment
% to a slack row or column.  I.e., when there are two or more missing
% image points or model points.
%
% A problem with this modified Sinkhorn algorithm is the following.
% Suppose all rows except the slack row are normalized.  It is possible that
% a nonslack value which was previously maximum in its row and column to now
% have a value that is less than the slack value for that column. (This
% nonslack value will still be greater than the slack value for that
% row.)  The same sort of thing can happen when columns are normalized.
% Intuitivitly, this seems like a bad thing: nonslack values that start
% off as maximum in their row and column should remain maximum in their
% row and column throughout this iteration.  The current algorithm
% attempts to prevent this from happening as follows.  After performing
% row normalizations, the values in the slack row are set so that their
% ratio to the nonslack value in that column which was previously maximum
% is the same as this ratio was prior to row normalization.  A similar
% thing is done after column normalizations.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function normalizedMat = sinkhornImp(M)

iMaxIterSinkhorn=60;          % In PAMI paper
fEpsilon2 = 0.001;            % Used in termination of Sinkhorn Loop.

iNumSinkIter = 0;
[nbRows, nbCols] = size(M);

fMdiffSum = fEpsilon2 + 1;    % Set "difference" from last M matrix above
                              % the loop termination threshold

% Get the positions and ratios to slack of the nonslack elements that
% are maximal in their row and column.
[posmax, ratios] = maxPosRatio(M);



while(abs(fMdiffSum) > fEpsilon2 & iNumSinkIter < iMaxIterSinkhorn)
    Mprev = M;  % Save M from previous iteration to test for loop termination

    % Col normalization (except outlier row - do not normalize col slacks
    % against each other)
    McolSums = sum(M, 1);  % Row vector.
    McolSums(nbCols) = 1;  % Don't normalize slack col terms against each other.
    McolSumsRep = ones(nbRows,1) * McolSums ;
    M = M ./ McolSumsRep;

    % Fix values in the slack column.
    for i = 1:size(posmax,1)
      M(posmax(i,1),nbCols) = ratios(i,1)*M(posmax(i,1),posmax(i,2));
    end

    % Row normalization (except outlier row - do not normalize col slacks
    % against each other)
    MrowSums = sum(M, 2);  % Column vector.
    MrowSums(nbRows) = 1;  % Don't normalize slack row terms against each other.
    MrowSumsRep = MrowSums * ones(1, nbCols);
    M = M ./ MrowSumsRep;

    % Fix values in the slack row.
    for i = 1:size(posmax,1)
      M(nbRows,posmax(i,2)) = ratios(i,2)*M(posmax(i,1),posmax(i,2));
    end

    iNumSinkIter=iNumSinkIter+1;
    fMdiffSum=sum(abs(M(:)-Mprev(:)));

end

normalizedMat = M;

return;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% maxPosRatio    Get the positions and ratios of the maximum nonslack values.
%     [POS, RATIOS] = maxPosRatio(ASSIGNMAT)
%         ASSIGNMAT is an image to model assignment matrix.  The last
%         row and column of ASSIGNMAT are slack, and are used to indicate
%         no match for that row or column. POS returns an Mx2 matrix of the
%         positions in ASSIGNMAT that are maximal within the respective
%         row and column. M is the number of elements in ASSIGNMAT that are
%         maximum in a row and column, and may equal zero.  The Kth row
%         of POS gives [ROWPOS COLPOS], the row and column position of the
%         Kth maximal element.  RATIOS returns an Mx2 matrix of the ratios
%         of these maximal values to the slack values in those rows and
%         columns.  The Kth row of RATIOS is [RRATIO, CRATIO] where RRATIO
%         (CRATIO) is the respective row (column) slack value for the Kth
%         maximal element divided by the value of the Kth maximal element.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pos, ratios] = maxPosRatio(assignMat)

pos = [];
ratios = [];

nrows = size(assignMat,1);
ncols = size(assignMat,2);
nimgpnts  = nrows - 1;
nmodpnts = ncols - 1;

% Iterate over all columns of assignMat.
for k = 1 : nmodpnts
    [vmax imax] = max(assignMat(:,k));                  % Max value in column k.

    if imax == nrows
        continue;                       % Slack value is maximum in this column.
    end;

    % Check if the max value in the column is maximum within its row.
    if all(vmax > assignMat(imax,[1:k-1,k+1:ncols]))
        pos = [pos; [imax, k]];     % This value is maximal in its row & column.

        % Compute the ratios to row and column slack values.
        rr = assignMat(imax,ncols)/assignMat(imax,k);
        cr = assignMat(nrows,k)/assignMat(imax,k);
        ratios = [ratios; [rr cr]];
    end
end

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NUMMATCHES    How many image points have found a model point match?
%     NUM = numMatches(ASSIGNMAT)
%         ASSIGNMAT is an image to model assignment matrix.  The last
%         row and column of ASSIGNMAT are slack, and are used to indicate
%         no match for that row or column. Image point number I matches
%         model point number M if ASSIGNMAT(I,M) is strictly greater
%         than all other entries in row I and column M of ASSIGNMAT.
%         NUM returns the number of matching points.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function num = numMatches(assignMat)

num = 0;

nrows = size(assignMat,1);
ncols = size(assignMat,2);
nimgpnts  = nrows - 1;
nmodpnts = ncols - 1;

for k = 1 : nmodpnts
    [vmax imax] = max(assignMat(:,k));                  % Max value in column k.

    if imax == nrows
        continue;                       % Slack value is maximum in this column.
    end;

    if all(vmax > assignMat(imax,[1:k-1,k+1:ncols]))
        num = num + 1;              % This value is maximal in its row & column.
    end
end

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
