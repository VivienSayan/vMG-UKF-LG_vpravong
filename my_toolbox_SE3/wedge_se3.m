function XI = wedge_se3(xi)
% Input:
% xi = element of R^d (the angular variable(s) must be first)
%
% Output: 
% compute the isomorphism from R^d to the Lie algebra
%
% This operator is overloaded : the result depends on the length d

if size(xi,1) == 3
    XI = wedge_so3(xi);    
elseif size(xi,1) == 6
    XI = [wedge_so3(xi(1:3)) xi(4:6);...
          0       0         0     0 ];    
end

end