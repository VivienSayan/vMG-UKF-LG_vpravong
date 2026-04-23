function Y = h_lidar(chi,landmarks,index_visible_lm,r)
% Y = h_lidar(chi,landmarks,index_visible_lm,r)
% "Unscented Kalman Filtering on Lie Groups (M.Brossard, S.Bonnabel, J-P Condomines, IROS 2017)"

Y = [];
for i = index_visible_lm
    Yi = inv(chi)*[landmarks(:,i);1];
    Y = [Y;Yi(1:3)];
end
Y = Y + r;

end