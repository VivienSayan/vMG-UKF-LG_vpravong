function Ad_chi = Adjoint_SE3(chi)
% Input: 
% chi (element of the group SE(3))
%
% Output: compute the Adjoint such that
% exp( [Adjoint_SE3(chi) * xi]^ ) = chi * exp(xi^) * chi^(-1)
%
% ^ denotes the "wedge" operator -> see the function wedge_se3()

R = chi(1:3,1:3);
t = chi(1:3,4);
T = wedge_so3(t);

Ad_chi = [R    zeros(3,3);...
          T*R     R];


%% Unit test

% chi = [eul2rotm([1,2,3]), [3;2;1]; 0 0 0, 1];
% xi = [0.01; 0.02; 0.03; 0.01; 0.02; 0.03];

% expm( wedge_se3(Adjoint_SE3(chi)*xi) )
% chi*expm(wedge_se3(xi))*chi^(-1)

