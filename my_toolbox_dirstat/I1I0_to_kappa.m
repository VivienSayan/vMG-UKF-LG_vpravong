function kappa = I1I0_to_kappa(I1I0)
% kappa = I1I0_to_kappa(I1I0)
% Using the method of "Modeling Data using Directional Distributions, Dhillon and Sra, 2003, equation (3.19)"
% Compute kappa = (I1I0*2-I1I0^3)/(1-I1I0^2)
%
% Inputs:
% I1I0: ratio of Bessel functions I_1(kappa)/I_0(kappa)
%
% Outputs:
% kappa: concentration parameter associated with the ratio of Bessels

kappa = (I1I0*2-I1I0^3)/(1-I1I0^2);

end