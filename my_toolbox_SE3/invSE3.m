function inv_chi = invSE3(chi)
% inv_chi = invSE3(chi)
% Input:
% chi: element of SE(3)
% 
% Output:
% inv_chi: inverse of chi
inv_chi = eye(size(chi));
R = chi(1:3,1:3);
inv_R = R';
inv_chi(1:3,1:3) = inv_R;
inv_chi(1:3,4:end) = -inv_R*chi(1:3,4:end);
end

