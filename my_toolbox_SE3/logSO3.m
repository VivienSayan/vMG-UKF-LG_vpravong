function omega = logSO3(R)
% res = logSO3(R)
% Input
% R : element of SO(3)
%
% Output:
% omega = [th1, th2, th3]
%
% OMEGA = [omega]^ = [0   -th3  th2]
%                    [th3  0   -th1]
%                    [-th2 th1  0  ]

theta = acos((trace(R)-1)/2) ;
if theta == 0
    OMEGA = zeros(3);
else
    OMEGA = theta/(2*sin(theta))*(R-R');
end
omega = [-OMEGA(2,3);OMEGA(1,3);-OMEGA(1,2)];
end
