function error = computeErrorSE3(traj,trajRef,i)
% error = computeErrorSE3(traj,trajRef,i)
%
% Inputs:
% traj : logs estimated trajectory
% traRef : logs reference trajectory
% i : last iteration
%
% Outputs:
% error : contains different error measures
% error.errorRot (deg)
% error.errorPos (m)
% error.errorLogPos (m)
% error.errorX (m)
% error.errorY (m)
% error.errorZ (m)

errorRot = zeros(i,1);
errorPos = zeros(i,1);
errorLogPos = zeros(i,1);
errorX = zeros(i,1);
errorY = zeros(i,1);
errorZ = zeros(i,1);

for j = 1:i
    RotRef = trajRef.Rot(:,:,j); Rot = traj.Rot(:,:,j);
    errorRot(j) = norm(logSO3(RotRef'*Rot));

    chi_error_j = invSE3(trajRef.chi(:,:,j))*traj.chi(:,:,j);
    tmp = logSE3(chi_error_j);
    errorLogPos(j) = norm(tmp(4:6));
    errorPos(j) = norm(traj.trans(:,j)-trajRef.trans(:,j));
    errorX(j) = norm(traj.trans(1,j)-trajRef.trans(1,j));
    errorY(j) = norm(traj.trans(2,j)-trajRef.trans(2,j));
    errorZ(j) = norm(traj.trans(3,j)-trajRef.trans(3,j));
end

error.errorRot = errorRot *180/pi;
error.errorPos = errorPos;
error.errorLogPos = errorLogPos;
error.errorX = errorX;
error.errorY = errorY;
error.errorZ = errorZ;
end