function traj = initialTrajSE3(eul_angles0,trans0,NbSteps,varargin)
% traj = initialTrajSE3(eul_angles0,trans0,NbSteps,varargin)
% Inputs:
% eul_angles0: initial Euler angles (in rad) -> [Yaw (around ez), Pitch (around ey), Roll (around ex)]
% trans0: initial translation vector (X,Y,Z)
% NbSteps: allocate number of timesteps 
% varargin: additional translational vectors (optional)
%
% Output:
% traj: initialization of the trajectory
% traj.eul_angles
% traj.Rot
% traj.trans
% traj.chi

traj.eul_angles = zeros(3,NbSteps); traj.eul_angles(:,1) = eul_angles0;
traj.Rot = zeros(3,3,NbSteps); traj.Rot(:,:,1) = eul2rotm(eul_angles0');
traj.trans = zeros(3,NbSteps); traj.trans(:,1) = trans0;

if (nargin == 4)
    add_pos0 = varargin{1};
    Nb_pos = size(add_pos0,2);
    traj.add_pos = zeros(3,Nb_pos,NbSteps); traj.add_pos(:,:,1) = add_pos0;

    traj.chi = zeros(4+Nb_pos,4+Nb_pos,NbSteps);
    traj.chi(:,:,1) = state2chiSE3(eul_angles0,trans0,add_pos0);
else
    traj.chi = zeros(4,4,NbSteps);
    traj.chi(:,:,1) = state2chiSE3(eul_angles0,trans0);
end