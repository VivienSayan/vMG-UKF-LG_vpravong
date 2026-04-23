function xi = logSE3(chi)
% xi = logSE3(chi) (See Barfoot & Furgale, Associating uncertainty for SE(3)... (2014), appendix B)
%
% Input:
% chi : element of SE(3)
% 
% Output: 
% xi : element of R^6
% N = number of additional translational vectors
% xi = [xi_angles, xi_trans];

C = chi(1:3,1:3);
r = chi(1:3,4);
[vecs,eigs] = eig(C);
[~,idx] = min(abs(diag(eigs)-1));
a = vecs(:,idx);

phi = acos(1/2*(trace(C)-1));
if phi == 0
    a = zeros(3,1);
    iJ = eye(3);
else
    A = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
    a_ta = a*transpose(a);
    C1 = cos(phi)*eye(3) + (1-cos(phi))*a_ta + sin(phi)*A;
    C2 = cos(-phi)*eye(3) + (1-cos(-phi))*a_ta + sin(-phi)*A;
    if norm(C2-C)<norm(C1-C)
        phi = -phi;
    end
    iJ = phi/2*cot(phi/2)*eye(3) + (1-phi/2*cot(phi/2))*a_ta - phi/2*A;
end
rho = iJ*r;
xi = real([phi*a;rho]);
end
