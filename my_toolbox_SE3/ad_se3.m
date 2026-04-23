function ad_XI = ad_se3(XI)
% See "Barfoot, 2014, associating uncertainty..." (eq. 12)
% Input (element of the Lie algebra se(3)):
% XI =  [0   -th3  th2  xi_x]
%       [th3  0   -th1  xi_y]
%       [-th2 th1  0    xi_z]
%       [0    0    0     0  ]
%
% Output: compute the "small" adjoint such that
% [ad_se3(XI_a) * xi_b]^ = xi_a^ * xi_b^ - xi_b^ * xi_a^
%
% ^ denotes the "wedge" operator -> see the function wedge_se3()

% xi = [th1, th2, th3, xi_x, xi_y, xi_z]
xi = [-XI(2,3); XI(1,3); -XI(1,2); XI(1,4); XI(2,4); XI(3,4)]; 
phi = xi(1:3);
rho = xi(4:6);
XI_phi = wedge_so3(phi);
XI_rho = wedge_so3(rho);

% Since xi first contains the angular components, then the skew of the rho
% part is located on the down-left side
ad_XI = [XI_phi  zeros(3,3);...
         XI_rho  XI_phi];
end

%% Unit test

% xi_a = [0.01; 0.02; 0.03; 0.01; 0.02; 0.03];
% xi_b = cos(xi_a);

% wedge_se3( ad_se3(wedge_se3(xi_a))*xi_b )
% wedge_se3(xi_a)*wedge_se3(xi_b) - wedge_se3(xi_b)*wedge_se3(xi_a)