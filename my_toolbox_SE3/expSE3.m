function chi = expSE3(xi)
% chi = expSE3(xi)
% Input:
% xi: element of R^6
%
% if (ex,ey,ez) are the body-fixed unit vectors
% xi = [xi_th1 (around ex), xi_th2 (around ey), xi_th3 (around ez) ,xi_x, xi_y, xi_z]
% 
% Output:
% chi: element of the group SE(3)

angles_norm = norm(xi(1:3));
if(angles_norm == 0)
    chi = [eye(3)      xi(4:6);...
           zeros(1,3)     1] ;
else
    XI = [wedge_so3(xi(1:3))    xi(4:6);... % in Lie algebra
             0     0     0         0];
    chi = eye(4) + XI + 1/angles_norm^2*(1-cos(angles_norm))*XI^2 + 1/angles_norm^3*(angles_norm-sin(angles_norm))*XI^3; % in Lie group
end
end

%% Unit test

% xi = [0.01; 0.02; 0.03; 0.5; 0.2; 0.32];

% expm(wedge_se3(xi))
% expSE3(xi)