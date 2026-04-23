function M = moyenne_diag_k(MC_Covs)

    % Dimensions
    [d1,d2,kmax,Ns] = size(MC_Covs);
    if d1 ~= 6 || d2 ~= 6
        error('MC_Covs doit être de taille 6 x 6 x kmax x Ns');
    end

    % Initialisation
    M = zeros(6, kmax);

    % Boucles explicites
    for k = 1:kmax
        for i = 1:6
            M(i,k) = mean( MC_Covs(i,i,k,:) );
        end
    end
end
