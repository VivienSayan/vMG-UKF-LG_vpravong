function [ym,index_visible_lm,dimy,sqrtR,sqrtR_true] = simul_meas_hlidar(chi,landmarks,R_true,R,rstream)
% [ym,index_visible_lm,dimy,sqrtR] = simul_meas_hlidar(chi,landmarks,R_true,R,rstream)
% Inputs:
% chi : element of SE(3)
% landmarks: table of size 3xm (for m landmarks)
% R_true: true covariance matrix of the noise
% R : estimated covariance matrix of the noise
% rstream : RNG stream
%
% Outputs:
% Simulate cartesian coordinates [x_i,y_i,z_i] of the landmark i with respect to the robot's body frame
% ym : observed measurement
% index_visible_lm: index of the observed landmarks
% dimy: dimension of the measurement vector
% sqrtR: cholesky factor of R (augmented dimension)
% sqrtR_true: cholesky factor R_true (augmented dimension)

Rot = chi(1:3,1:3);
X = chi(1:3,4);

sqrtR_true = chol(R_true,'upper');
r_x = R(1,1); r_y = R(2,2); r_z = R(3,3);

ym = []; diag_sqrtR = []; index_visible_lm = []; diag_sqrtR_true = [];
for ind = 1:size(landmarks,2)
    lm = landmarks(:,ind); % landmarks pos in inertial frame
    dist = norm(X-lm);
    if (0.01 <= dist && dist <= 120)
        yi = invSE3(chi)*[lm;1]; % convert lm pos from inertial to body frame
        yi = yi(1:3)+sqrtR_true*randn(rstream,3,1);
        ym = [ym;yi];
        diag_sqrtR = [diag_sqrtR sqrt(r_x) sqrt(r_y) sqrt(r_z)];
        diag_sqrtR_true = [diag_sqrtR_true sqrt(R_true(1,1)) sqrt(R_true(2,2)) sqrt(R_true(3,3))];
        index_visible_lm = [index_visible_lm ind];
    end
end

dimy = size(ym,1);
sqrtR = diag(diag_sqrtR);
sqrtR_true = diag(diag_sqrtR_true);

end