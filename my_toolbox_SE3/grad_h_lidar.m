function H = grad_h_lidar(attitude,position,landmarks,index_visible_lm)
% H = grad_h_lidar(attitude,position,landmarks,index_visible_lm)
% Inputs:
% attitude: rotation matrix in SO(3)
% position: position vector in R^3
% landmarks: table of size 3xm (for m landmarks)
% dim(index_visible_lm): table of size 1x max_nb
%
% Output:
% H: gradient of hlidar observation model

R = attitude; x = position;
E1 = [0 0 0; 0 0 -1; 0 1 0]; E2 = [0 0 1; 0 0 0; -1 0 0]; E3 = [0 -1 0; 1 0 0; 0 0 0];

H = [];
for i = 1:size(index_visible_lm,2)
    bi = landmarks(:,i);
    Hi = [R'*E1*(bi-x), R'*E2*(bi-x), R'*E3*(bi-x), -R'] ;
    H = [H;Hi];
end


end