function chi = exp_multiSE3(xi)
% chi = exp_multiSE3(xi)
% Input:
% xi: element of R^{6+3*N}
%
% if (ex,ey,ez) are the body-fixed unit vectors
% xi = [xi_th1 (around ex), xi_th2 (around ey), xi_th3 (around ez) ,xi_x, xi_y, xi_z, xi_x_1,xi_y_1,xi_z_1, ... , xi_x_N, xi_y_N, xi_z_N]
% 
% Output:
% chi: element of the group SE_{1+N}(3)
% N = number of additional translational vectors

xi_angles = xi(1:3);
xi_trans = xi(4:end);
angles_norm = norm(xi_angles);
Nplus1 = length(xi)/3-1; % total number of translational vectors (body translational vector + additional translational vectors)
if(angles_norm == 0)
    chi = eye(3+Nplus1);
    chi(1:3,4:end) = reshape(xi_trans,[3 Nplus1]);
else
    XI = zeros(3+Nplus1);
    XI(1:3,1:3) = wedge_so3(xi_angles);
    XI(1:3,4:end) = reshape(xi_trans,[3 Nplus1]);
    chi = eye(3+Nplus1) + XI + 1/angles_norm^2*(1-cos(angles_norm))*XI^2 + 1/angles_norm^3*(angles_norm-sin(angles_norm))*XI^3;
end
end

%% Unit test

% xi = [0.01; 0.02; 0.03; 0.5; 0.2; 0.32];

% expm(wedge_se3(xi))
% exp_multiSE3(xi)