function J_SE3 = JacSE3(XI,filter)
% (See Barfoot,2014,"Associating uncertainty..." Appendix)
% compute the RIGHT jacobian if filter = LEFT
% compute the LEFT jacobian if filter = RIGHT
%
% Inputs:
% xi = [th1, th2, th3, ,xi_x, xi_y, xi_z]
% XI =  [0   -th3  th2  xi_x]   element of the Lie algebra
%       [th3  0   -th1  xi_y]
%       [-th2 th1  0    xi_z]
%       [0    0    0     0  ]
% filter: LEFT or RIGHT UKF
%
% Output:
% J_SE3 = group jacobian

xi = [-XI(2,3); XI(1,3); -XI(1,2); XI(1,4); XI(2,4); XI(3,4)];

%Compute Right Jacobian
taille_xi = length(xi);
tolerance = 1e-12;
Nb_add_pos = taille_xi/3-2;
phi = xi(1:3);

ph = norm(phi);
if ph < tolerance % If the angle is small, fall back on the series representation
    switch filter
        case "LEFT"
            J_SE3 = JacSE3_Right_Series(XI,10);
        case "RIGHT"
            J_SE3 = JacSE3_Left_Series(XI,10);
        otherwise
            error("filter does not exist")
    end
    % Q = xi2Q_Series(xi(1:6),10,filter);
    % J_SE3 = [J_SO3   zeros(3,3);...
    %            Q     J_SO3];
else
    axis = phi/ph;
    cph = (1 - cos(ph))/ph;
    sph = sin(ph)/ph;
    J_SO3 = sph * eye(3) + (1 - sph) * (axis)*(axis)' + cph * wedge_so3(axis);

    Q = xi2Q(xi,filter);
    if Nb_add_pos > 0
        J_SE3 = kron(eye(Nb_add_pos+2),J_SO3);
        J_SE3(1:3,4:6) = Q;
        for i = 1:Nb_add_pos
            Q = xi2Q([phi;xi(4+3*i:6+3*i)],filter);
            J_SE3(1:3,4+3*i:6+3*i) = Q; % :/
        end
    else
        J_SE3 = [J_SO3  zeros(3,3);...
                 Q      J_SO3];
    end

end

% must respect properties:
%det(J_SO3)^(2+Nb_add_pos) =
%det(J_SE3)

end