function xi = vee_se3(XI)
% xi = vee_se3(XI)
% Input: 
% XI : element of se(3)
% XI =  [0   -th3  th2  xi_x]
%       [th3  0   -th1  xi_y]
%       [-th2 th1  0    xi_z]
%       [0    0    0     0  ]
%
% Output: compute the isomorphism from se(3) to R^6
% xi : element of R^6

xi = [-XI(2,3); XI(1,3); -XI(1,2); XI(1,4); XI(2,4); XI(3,4)];
end