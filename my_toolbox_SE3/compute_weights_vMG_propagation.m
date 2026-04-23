function [Wm,Wc,ksi_x,ksi_th1,ksi_th2,ksi_th3,KSI] = compute_weights_vMG_propagation(n,m,ut_alpha,ut_kappa,kappa_th1,kappa_th2,kappa_th3,kappa_q_th1,kappa_q_th2,kappa_q_th3)
% compute horwood weights and nodes based on UT parameters
% ut_alpha,ut_kappa, vm_kappa_phi1, vm_kappa_phi2 and vm_kappa_phi3,
% vm_kappa_q1,vm_kappa_q2,vm_kappa_q3
% The state belongs to R^n x S^m.

%-- translational variables
lambda = ut_alpha^2*(n+ut_kappa)-n;
ksi_x = sqrt(n+lambda);
W_x = 1/(2*(n+lambda));

%-- rotational variables
% states khi
B1 = 1-bessel_ratio(0,kappa_th1,10); 
B2 = 1-bessel_ratio(0,kappa_th1,10)*bessel_ratio(1,kappa_th1,10);
ksi_th1 = acos(B2/(2*B1)-1);
W_th1 = B1^2/((4*B1-B2));

B1 = 1-bessel_ratio(0,kappa_th2,10); 
B2 = 1-bessel_ratio(0,kappa_th2,10)*bessel_ratio(1,kappa_th2,10);
ksi_th2 = acos(B2/(2*B1)-1);
W_th2 = B1^2/((4*B1-B2));

B1 = 1-bessel_ratio(0,kappa_th3,10); 
B2 = 1-bessel_ratio(0,kappa_th3,10)*bessel_ratio(1,kappa_th3,10);
ksi_th3 = acos(B2/(2*B1)-1);
W_th3 = B1^2/((4*B1-B2));

% input noises
B1 = 1-bessel_ratio(0,kappa_q_th1,10); 
B2 = 1-bessel_ratio(0,kappa_q_th1,10)*bessel_ratio(1,kappa_q_th1,10);
ksi_q_th1 = acos(B2/(2*B1)-1);
W_q_th1 = B1^2/((4*B1-B2));

B1 = 1-bessel_ratio(0,kappa_q_th2,10); 
B2 = 1-bessel_ratio(0,kappa_q_th2,10)*bessel_ratio(1,kappa_q_th2,10);
ksi_q_th2 = acos(B2/(2*B1)-1);
W_q_th2 = B1^2/((4*B1-B2));

B1 = 1-bessel_ratio(0,kappa_q_th3,10); 
B2 = 1-bessel_ratio(0,kappa_q_th3,10)*bessel_ratio(1,kappa_q_th3,10);
ksi_q_th3 = acos(B2/(2*B1)-1);
W_q_th3 = B1^2/((4*B1-B2));
%--

W0 = 1 -2*W_th1-2*W_th2-2*W_th3 -2*n*W_x -2*W_q_th1-2*W_q_th2-2*W_q_th3;

Wm = [W0 W_th1 W_th2 W_th3 repmat(W_x,1,n) W_q_th1 W_q_th2 W_q_th3 W_th1 W_th2 W_th3 repmat(W_x,1,n) W_q_th1 W_q_th2 W_q_th3];
Wc = Wm;

KSI_X = ksi_x*eye(3);
KSI_TH = diag([ksi_th1;ksi_th2;ksi_th3]);
KSI_Q_TH = diag([ksi_q_th1;ksi_q_th2;ksi_q_th3]);
KSI = blkdiag(KSI_TH,KSI_X,KSI_Q_TH,KSI_X);

end