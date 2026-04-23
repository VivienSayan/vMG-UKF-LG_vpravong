function [Rot,eul_angles,trans,add_trans] = chi2stateSE3(chi)
% Input:
% chi : element of the group SE(3)
%
% Output:
% Rot : rotation matrix
% eul_angles : Euler angles associated with the rotation matrix
% trans : translational vector
% add_pos : additional translational vectors
%
% [Rot, eul_angles, trans, add_pos] = chi2stateSE3(chi)

Rot = chi(1:3,1:3);
eul_angles = rotm2eul(Rot)'; % rotm2eul returns 1x3 array, so transpose it to get 3x1 array
trans = chi(1:3,4);
if size(chi,2) > 4
    add_trans = chi(1:3,5:end);
else
    add_trans = [];
end
end
