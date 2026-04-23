function xi = vee_so3(XI)
% xi = vee_so3(XI)
% XI belongs to so(3) and xi belongs to R^3
% XI =  [0   -th3  th2 ]
%       [th3  0   -th1 ]
%       [-th2 th1  0   ]
xi = [-XI(2,3); XI(1,3); -XI(1,2)];
end