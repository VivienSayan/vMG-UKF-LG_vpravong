function chi = state2chiSE3(eul_angles,trans,varargin)
% chi = state2chiSE3(eul_angles,trans,varargin)
% arg1 = eul_angles (in rad) -> Yaw (around Z), Pitch (around Y), Roll (around X)
%        Rot = eul2rotm(eul_angles) = Rz * Ry * Rx
% arg2 = trans (X,Y,Z)
% arg3 = varargin (additional poses, optional)

Rot = eul2rotm(eul_angles');
if(nargin == 3)
    add_pos = varargin{1};
    Nb_pos = size(add_pos,2);
    chi = eye(4+Nb_pos);
    chi(1:3,1:4) = [Rot trans];
    if Nb_pos > 0
        chi(1:3,5:end) = add_pos;
    end
else
    chi = [Rot trans; zeros(1,3) 1];
end

end
