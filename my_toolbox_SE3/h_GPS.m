function y = h_GPS(chi,r)
% y = h_GPS(chi,r)
% See "Unscented Kalman Filtering on Lie Groups (M.Brossard, S.Bonnabel, J-P Condomines, IROS 2017)"

y = chi*[0;0;0;1] + [r;0];
y = y(1:3);

end