function y = bessel_ratio(nu,kappa,N)
% y = bessel_ratio(nu,kappa,N)
% Compute I_{nu+1}(kappa) / I_{nu}(kappa) using the method of 
% "Recursive Nonlinear Filtering for Angular Data  Based on Circular
% Distributions" Kurz, Hanebeck, 2013
%
% Inputs:
% nu: order of the Bessel function
% kappa: concentration parameter
% N: number of iterations
%
% Output:
% y: ratio of Bessel functions

o = min([nu,10]);
r = nan(1,N+1);
for i = 0:N
    r(i+1) = kappa/( o+i+0.5+ sqrt((o+i+1.5)^2+kappa^2) );
end
for i = 1:N
    for k = 0:N-1
        r(k+1) = kappa/( o+k+1+ sqrt((o+k+1)^2+kappa^2*r(k+2)/r(k+1)) );
    end
end
y = r(1);
i = o;
while i > nu
    y = 1/(2*i/kappa+y);
    i = i-1;
end
end