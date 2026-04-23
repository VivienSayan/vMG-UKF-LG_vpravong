function J = JacSE3_Right_Series( XI, N )
% J = JacSE3_Right_Series( XI, N )
% (This function is overloaded: J is wether the right jacobian for SE(3) or SO(3) depending on the size of XI)
% 
% Input
% XI : element of the Lie algebra
% N: order of the serie
%
% Output:
% J : right jacobian by series expansion

if size(XI,1) == 3 
    
    J = eye(3);
    adn = eye(3);
    ad = -XI; % in general, it should be "little adjoint" operator but for so(3), adjoint(xi^) = xi^
    for n = 1:N
        adn = adn*ad/(n + 1);    
        J = J + adn;
    end
    
elseif size(XI,1) == 4    
    
    J = eye(6);
    adn = eye(6);
    ad = -ad_se3(XI);
    for n = 1:N
        adn = adn*ad /(n + 1);    
        J = J + adn;
    end      
end
end