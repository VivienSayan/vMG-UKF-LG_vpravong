function XI = wedge_so3(xi)
% Input:
% xi = element of R^3 -> xi = [th1; th2; th3]
%
% Output: 
% compute the isomorphism from R^3 to the Lie algebra
%
% XI =  [0   -th3  th2]
%       [th3  0   -th1]
%       [-th2 th1  0  ]

XI = [0       -xi(3)   xi(2);...
      xi(3)     0     -xi(1);...
     -xi(2)   xi(1)     0];
end