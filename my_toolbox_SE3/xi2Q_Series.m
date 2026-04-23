function Q = xi2Q_Series(xi,N,filter)
% See Barfoot et al. (2014) Appendix (eq.102) 
switch filter
    case "RIGHT"
        phi = xi(1:3); % rotation part
        rho = xi(4:6); % translation part
    case "LEFT"
        phi = -xi(1:3); % rotation part
        rho = -xi(4:6); % translation part
    otherwise
        error("filter does not exist")
end

sum2 = 0;
for n=0:N
    sum1 = 1/factorial(n+2) * wedge_so3(phi)^n * wedge_so3(rho);
    for m=1:N
        sum1 = sum1 + 1/factorial(n+m+2) * wedge_so3(phi)^n * wedge_so3(rho) * wedge_so3(phi)^m;
    end
    sum2 = sum2 + sum1;
end

Q = sum2;

end