function traj = updateTrajSE3(traj,chi,Rot,eul_angles,trans,k,varargin)
% traj = updateTrajSE3(traj,chi,Rot,eul_angles,trans,k,varargin)

traj.eul_angles(:,k) = eul_angles;
traj.trans(:,k) = trans;
traj.Rot(:,:,k) = Rot;
traj.chi(:,:,k) = chi;
if (nargin == 7)
    add_pos = varargin{1};
    traj.add_pos(:,:,k) = add_pos;
end

end