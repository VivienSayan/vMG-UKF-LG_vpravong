clear all;
%clc; 
%close all;

addpath my_toolbox_dirstat;
addpath my_toolbox_SE3;
addpath data_3DPose;

%data = 'data_3DPose/traj_eight_SE3.mat';
%data = 'data_3DPose/robot_rpi_meas2.mat';
data = 'data_3DPose/Euro_dataset_V1_03_difficult_velocities.mat'; 
load(data);
load('curveFit.mat');
%load('curveFit_power2_26_to_90.mat'); 
%var_to_kappa = var_to_kappa_power2_26_to_90;

% Select hyper parameter
ut_type = 'ut_gvm' % 'ut_sr' / 'ut_standard' / 'ut_gvm'
n = 3; % number of gaussian variables in the state
m = 3; % number of circular variables in the state
filter_side = 'LEFT' % 'LEFT' / 'RIGHT'
correl_prediction = 'y' % 'y'/'n' (take correlation Pxth into account or not in the prediction step)
correl_correction = 'y' % 'y'/'n' (take correlation Pxth into account or not in the correction step)
bool = 1; % 1 (apply update step) / nan (no update step)
CLVQ = 'no'; % 'yes' / 'no'
    gamma_th1 = 1/10; gamma_th2 = 1/10; gamma_th3 = 1/10; gamma_x = 1/20; gamma_y = 1/20; gamma_z = 1/20;
    N = 100; ndim = [1;2;3];
CLVQbis =  'no'; % 'yes' / 'no'
    gamma_th1_bis = 1/10; gamma_th2_bis = 1/10; gamma_th3_bis = 1/10; gamma_x_bis = 1/20; gamma_y_bis = 1/20; gamma_z_bis = 1/20;
    Nbis = 100; ndimbis = [1;2;3];
sensor = 'Lidar' % 'GPS' / 'Lidar' / 'full_state'
f_obs = 1 % measurement frequency (Hz)

% Monte-carlo simulations
seed = 23; % Choose seed
rstream = RandStream('dsfmt19937','Seed',seed);
Ns = 1; % Number of Monte-Carlo Simulations

% Landmarks position (expressed in inertial frame)
landmarks = [ 0  4 0;...
            -0.5 1 2;...
              0  0 0];

% measurement acquisition
obs = zeros(1,kmax);
dt_obs = 1/(dt*f_obs);
for i = 1:round(dt_obs):kmax
    obs(i) = 1;
end

E1 = [0 0 0; 0 0 -1; 0 1 0]; E2 = [0 0 1; 0 0 0; -1 0 0]; E3 = [0 -1 0; 1 0 0; 0 0 0];

% set noise level
sigma_model = linspace(1,70,10).^2;

% Initialize ground truth
trajGT = initialTrajSE3([trajGT.psi(1);trajGT.theta(1);trajGT.phi(1)],[trajGT.x(:,1)],kmax);
for level = 10:10
    level

    % Allocate data logs---------------
    MC_results = struct('eul_angles',cell(1,Ns),'trans',cell(1,Ns),'Rot',cell(1,Ns),'chi',cell(1,Ns));
    MC_GT = struct('eul_angles',cell(1,Ns),'trans',cell(1,Ns),'Rot',cell(1,Ns),'chi',cell(1,Ns));
    MC_Covs = zeros(6,6,kmax,Ns);
    MC_computation_time = zeros(1,Ns);
    MC_rmse_rot = zeros(1,Ns);
    MC_rmse_pos = zeros(1,Ns);
    MC_I1I0_hat = zeros(3,kmax,Ns);
    MC_kappa_hat = zeros(3,kmax,Ns);
    MC_LG_Error_squared_k = zeros(6,kmax,Ns);
    %---------------------------------

    for mc = 1:Ns
        % initial state (true)
        eul_angles_GT = trajGT.eul_angles(:,1);
        x_GT = trajGT.trans(:,1);
        chiGT = state2chiSE3(eul_angles_GT,x_GT);
        [Rot_GT,eul_angles_GT,x_GT,add_pos] = chi2stateSE3(chiGT);
        trajGT = updateTrajSE3(trajGT,chiGT,Rot_GT,eul_angles_GT,x_GT,1);

        % initial covariance
    p0Roll = (60*pi/180)^2; p0Pitch = (60*pi/180)^2; p0Yaw = (60*pi/180)^2; p0Rot = [p0Roll;p0Pitch;p0Yaw];
        kappa_th1_hat = var_to_kappa(p0Roll); kappa_th2_hat = var_to_kappa(p0Pitch); kappa_th3_hat = var_to_kappa(p0Yaw);
        I1I0_th1_hat = bessel_ratio(0,kappa_th1_hat,10); I1I0_th2_hat = bessel_ratio(0,kappa_th2_hat,10); I1I0_th3_hat = bessel_ratio(0,kappa_th3_hat,10);
    p0x = (0.5)^2; p0y = (0.5)^2; p0z = (0.5)^2; p0Pos = [p0x;p0y;p0z];
        P = diag([p0Rot;p0Pos]); dim_state = size(P,1); S = chol(P,'upper');
        % initial estimate
        roll0 = eul_angles_GT(3) + circ_vmrnd(0,kappa_th1_hat,1);%sqrt(p0Roll)*randn(rstream,1,1);% 
        pitch0 = eul_angles_GT(2) + circ_vmrnd(0,kappa_th2_hat,1);%sqrt(p0Pitch)*randn(rstream,1,1);% 
        yaw0 = eul_angles_GT(1) + circ_vmrnd(0,kappa_th3_hat,1);%sqrt(p0Yaw)*randn(rstream,1,1);% 
        eul_angles = [yaw0;pitch0;roll0];
        x = x_GT + sqrtm(diag(p0Pos))*randn(rstream,size(p0Pos,1),1);
        chi = state2chiSE3(eul_angles,x);

        % process noise
        % true noise
        q_th1_true = sigma_model(level)*(1*pi/180)^2; kappa_q_th1_true = var_to_kappa(q_th1_true);
        q_th2_true = sigma_model(level)*(1*pi/180)^2; kappa_q_th2_true = var_to_kappa(q_th2_true);
        q_th3_true = sigma_model(level)*(1*pi/180)^2; kappa_q_th3_true = var_to_kappa(q_th3_true);
        q_x_true = sigma_model(level)*(0.01)^2;
        q_y_true = sigma_model(level)*(0.01)^2;
        q_z_true = sigma_model(level)*(0.01)^2;
        Q_true = diag([q_th1_true;q_th2_true;q_th3_true;q_x_true;q_y_true;q_z_true]);
        % kalman parameters
        q_th1 = sigma_model(level)*(1*pi/180)^2; kappa_q_th1 =  var_to_kappa(q_th1);%kappa_q_th1_true;
        q_th2 = sigma_model(level)*(1*pi/180)^2; kappa_q_th2 =  var_to_kappa(q_th2);%kappa_q_th2_true;
        q_th3 = sigma_model(level)*(1*pi/180)^2; kappa_q_th3 =  var_to_kappa(q_th3);%kappa_q_th_true;
        q_x = sigma_model(level)*(0.01)^2;
        q_y = sigma_model(level)*(0.01)^2;
        q_z = sigma_model(level)*(0.01)^2;
        Q = diag([q_th1;q_th2;q_th3;q_x;q_y;q_z]); dimq = size(Q,1);
        sqrtQ = chol(Q,'upper');

        switch sensor
            case 'GPS'
                r_x_true = sigma_model(1)*(0.001)^2; r_y_true = r_x_true; r_z_true = r_x_true; R_true = diag([r_x_true;r_y_true;r_z_true]);
                r_x = sigma_model(1)*(0.001)^2; r_y = r_x; r_z = r_x; R = diag([r_x;r_y;r_z]); 
            case 'Lidar'
                r_x_true = sigma_model(1)*(0.05)^2; r_y_true = r_x_true; r_z_true = r_x_true; R_true = diag([r_x_true;r_y_true;r_z_true]);
                r_x = sigma_model(1)*(0.05)^2; r_y = r_x; r_z = r_x; R = diag([r_x;r_y;r_z]);  
            otherwise
                error('no sensor');
        end
        dimy = size(R,1);
        sqrtR_true = chol(R_true,'upper');
        sqrtR = chol(R,'upper');

        % Set UT hyper-parameters
        UT_alpha = 1;
        UT_beta = 0;
        UT_kappa = nan;

        % Allocate data logs------------------------------------------
        trajukfLG = initialTrajSE3(eul_angles,x,kmax);
        Covs = zeros(dim_state,dim_state,kmax); Covs(:,:,1) = S'*S;
        I1I0_hat_log = zeros(3,kmax); I1I0_hat_log(:,1) = [I1I0_th1_hat;I1I0_th2_hat;I1I0_th3_hat];
        kappa_hat_log = zeros(3,kmax); kappa_hat_log(:,1) = [kappa_th1_hat;kappa_th2_hat;kappa_th3_hat];
        MC_Covs(:,:,1,mc) = S'*S;
        MC_I1I0_hat(:,1,mc) = [I1I0_th1_hat;I1I0_th2_hat;I1I0_th3_hat];
        MC_kappa_hat(:,1,mc) = [kappa_th1_hat;kappa_th2_hat;kappa_th3_hat];
        LG_Error_k = zeros(6,kmax); LG_Error_k(:,1) = logSE3(invSE3(chiGT)*chi);
        %--------------------------------------------------------------

        switch sensor
            case 'Lidar'
                [ym,index_visible_lm,dimy,sqrtR,sqrtR_true_tmp] = simul_meas_hlidar(trajGT.chi(:,:,1),landmarks,R_true,R,rstream);
                H = (-1) * grad_h_lidar(Rot_GT,x_GT,landmarks,index_visible_lm);
                Hk = zeros(9,6,kmax); Hk(1:dimy,1:6,1) = H;
                Rk = zeros(9,9,kmax); Rk(1:dimy,1:dimy,1) = sqrtR_true_tmp'*sqrtR_true_tmp;
                dim_obs_k = zeros(1,kmax); dim_obs_k(1) = dimy;
            case 'GPS'
                H = eye(dimy,dimy);
                Hk = zeros(dimy,dimy,kmax); Hk(1:dimy,1:dimy,1) = H;
                Rk = zeros(dimy,dimy,kmax); Rk(1:dimy,1:dimy,1) = R;
                dim_obs_k = zeros(1,kmax); dim_obs_k(1) = dimy;
        end

        for k = 2:kmax % ____________begin filtering_______________________
            start_filter = tic;
            dt = time(k)-time(k-1);

            % velocity inputs polluted by noise
            omega_th1 = u(1,k) + circ_vmrnd(0,kappa_q_th1_true,1);%sqrt(q_th1_true)*randn(rstream,1,1);%
            omega_th2 = u(2,k) + circ_vmrnd(0,kappa_q_th2_true,1);%sqrt(q_th2_true)*randn(rstream,1,1);%
            omega_th3 = u(3,k) + circ_vmrnd(0,kappa_q_th3_true,1);%sqrt(q_th3_true)*randn(rstream,1,1);%
            ux = u(4,k)        + sqrt(q_x_true)*randn(rstream,1,1);  
            uy = u(5,k)        + sqrt(q_y_true)*randn(rstream,1,1);
            uz = u(6,k)        + sqrt(q_z_true)*randn(rstream,1,1);
            % propagate true pose
            chi_km1 = chiGT;
            chiGT = chiGT*exp_multiSE3([omega_th1;omega_th2;omega_th3;ux;uy;uz]*dt);
            [Rot_GT,eul_angles_GT,x_GT,add_pos] = chi2stateSE3(chiGT);
            trajGT = updateTrajSE3(trajGT,chiGT,Rot_GT,eul_angles_GT,x_GT,k);

            % generate sigma-points
            switch ut_type
                case 'ut_sr'
                    S_aug = blkdiag(S,sqrtQ); dim_state_aug = size(S_aug,1); % augment the state with the noises
                    x_aug = zeros(dim_state_aug,1); % augmented vector (in Lie algebra)
                    UT_kappa = 3-dim_state_aug;
                    [Wm,Wc,ksi,KSI,lambda] = compute_weights(dim_state_aug,UT_alpha,UT_beta,UT_kappa); % Weights
                    SigPts_01 = [zeros(dim_state_aug,1) -KSI*eye(dim_state_aug,dim_state_aug) KSI*eye(dim_state_aug,dim_state_aug)];
                    SigPts = zeros(dim_state_aug,2*dim_state_aug+1);
                    for i = 1:2*dim_state_aug+1
                        SigPts(:,i) = x_aug(:) + S_aug'*SigPts_01(:,i); % do the transpose if "lower" cholesky
                    end
                case 'ut_gvm'
                    S_aug = blkdiag(S,sqrtQ); dim_state_aug = size(S_aug,1); % augment the state with the noises
                    state_aug = zeros(dim_state_aug,1); % augmented state (in Lie algebra)
                    n_aug = dim_state_aug-2*m; % augmented gaussian variables
                    m_aug = dim_state_aug-n_aug; % augmented circular variables
                    UT_kappa = 3-n_aug;
                    concentration_th1 = var_to_kappa(P(1,1)); concentration_th2 = var_to_kappa(P(2,2)); concentration_th3 = var_to_kappa(P(3,3));
                    switch correl_prediction
                        case 'y'
                            Stmp = S;
                            M = P(4:dim_state,1:3)*P(1:3,1:3)^(-1);
                            Sxx_tilde = S_aug(4:dim_state,4:dim_state);
                            U = M*Stmp(1:3,1:3)'; % => U*U' = M*Pth*M';
                            %-------- replace Pxx = Pxx + MPthM' ----------
                            for j = 1:3
                                Sxx_tilde = cholupdate(Sxx_tilde,U(:,j),'+'); % --->>>> Equivalent to Sxx_tilde = chol(P(4:dim_state,4:dim_state) + M*P(1:3,1:3)*M','upper') but more stable
                            end
                            %----------------------------------------------
                            S_aug(4:dim_state,4:dim_state) = Sxx_tilde;
                        case 'n'
                            %M = P(4:dim_state,1:3)*P(1:3,1:3)^(-1);
                            %S_aug(4:dim_state,4:dim_state) = chol(P(4:dim_state,4:dim_state)+M*P(1:3,1:3)*M','upper');
                    end
                    [Wm,Wc,ksi_x,ksi_th1,ksi_th2,ksi_th3,KSI] = compute_weights_vMG_propagation(n_aug,m_aug,UT_alpha,UT_kappa,concentration_th1,concentration_th2,concentration_th3,kappa_q_th1,kappa_q_th2,kappa_q_th3); % Horwood Weights
                    SigPts_01 = [zeros(dim_state_aug,1) -KSI*eye(dim_state_aug,dim_state_aug) KSI*eye(dim_state_aug,dim_state_aug)];
                    SigPts = zeros(dim_state_aug,2*dim_state_aug+1);
                    for i = 1:2*dim_state_aug+1
                        SigPts(1:m,i)                     = SigPts_01(1:m,i);
                        SigPts(m+1:dim_state,i)           = S_aug(m+1:dim_state,m+1:dim_state)'*SigPts_01(m+1:dim_state,i); % do the transpose if "lower" cholesky
                        %--
                        SigPts(dim_state+1:dim_state+m,i) = SigPts_01(dim_state+1:dim_state+m,i);
                        SigPts(dim_state+m+1:end,i)       = S_aug(dim_state+m+1:end,dim_state+m+1:end)'*SigPts_01(dim_state+m+1:end,i);
                    end
                case 'ut_standard'
                    P_aug = blkdiag(P,Q); dim_state_aug = size(P_aug,1); 
                    state_aug = zeros(dim_state_aug,1);
                    UT_kappa = 3-dim_state_aug;
                    [Wm,Wc,ksi,KSI,lambda] = compute_weights(dim_state_aug,UT_alpha,UT_beta,UT_kappa);
                    SigPts_01 = [zeros(dim_state_aug,1) -KSI*eye(dim_state_aug,dim_state_aug) KSI*eye(dim_state_aug,dim_state_aug)];
                    SigPts = zeros(dim_state_aug,2*dim_state_aug+1);
                    for i = 1:2*dim_state_aug+1
                        SigPts(:,i) = state_aug(:) + sqrtm(P_aug)*SigPts_01(:,i);
                    end
                otherwise
                    error('UT type undefined')
            end

            switch CLVQ
                case 'yes'
                    gamma = diag([gamma_th1;gamma_th2;gamma_th3;gamma_x;gamma_y;gamma_z])*S; P = S'*S;
                    [SigPts,w_QO] = QO(gamma(ndim,ndim),N,x_aug,P,SigPts,ndim);
                case 'no'
                    ;
                otherwise
                    error("CLVQ -> you have to specify 'yes' or 'no' ")
            end

            chi_previous = chi; % save anterior pose
            %-------------------------- PREDICTION ----------------------------        
            % mean propagation without noise
            vel = [u(1,k);u(2,k);u(3,k);u(4,k);u(5,k);u(6,k)];
            chi = chi*exp_multiSE3(vel*dt); % -> propagated mean state in the Lie group
            chi_inv = invSE3(chi);

            for j = 2:2*dim_state_aug+1
                ksi_j = SigPts(1:dim_state,j);
                qj = SigPts(dim_state+1:end,j);
                uj = vel;
                switch filter_side
                    case 'LEFT'
                        chi_j = chi_previous*exp_multiSE3(ksi_j);
                        chi_j = chi_j*exp_multiSE3((uj+qj)*dt);
                        Xi_j = chi_inv*chi_j;
                    case 'RIGHT'
                        chi_j = exp_multiSE3(ksi_j)*chi_previous;
                        chi_j = chi_j*exp_multiSE3((uj+qj)*dt);
                        Xi_j = chi_j*chi_inv;
                    otherwise
                        error('filter does not exist');
                end
                SigPts(1:dim_state,j) = log_multiSE3(Xi_j); % propagated sigma-point j in the Lie algebra
            end

            % compute a priori estimate and covariance
            switch ut_type
                case 'ut_sr'
                    WSigPts = sqrt(Wc(2:end)).*SigPts(1:dim_state,2:end);
                    [~,RSx] = qr(WSigPts');
                    S = RSx(1:dim_state,1:dim_state);
                    P = S'*S;
                    I1I0_th1_hat = sqrt( sum(Wc.*cos(SigPts(1,:)),2)^2 + sum(Wc.*sin(SigPts(1,:)),2)^2 );
                    kappa_th1_hat = I1I0_to_kappa(I1I0_th1_hat);
                    I1I0_th2_hat = sqrt( sum(Wc.*cos(SigPts(2,:)),2)^2 + sum(Wc.*sin(SigPts(2,:)),2)^2 );
                    kappa_th2_hat = I1I0_to_kappa(I1I0_th2_hat);
                    I1I0_th3_hat = sqrt( sum(Wc.*cos(SigPts(3,:)),2)^2 + sum(Wc.*sin(SigPts(3,:)),2)^2 );
                    kappa_th3_hat = I1I0_to_kappa(I1I0_th3_hat);
                case 'ut_gvm'
                    WSigPts = sqrt(Wc(2:end)).*SigPts(1:dim_state,2:end);
                    [~,RSx] = qr(WSigPts');
                    S = RSx(1:dim_state,1:dim_state);
                    P = S'*S;
                    I1I0_th1_hat = sqrt( sum(Wc.*cos(SigPts(1,:)),2)^2 + sum(Wc.*sin(SigPts(1,:)),2)^2 );
                    kappa_th1_hat = I1I0_to_kappa(I1I0_th1_hat);
                    I1I0_th2_hat = sqrt( sum(Wc.*cos(SigPts(2,:)),2)^2 + sum(Wc.*sin(SigPts(2,:)),2)^2 );
                    kappa_th2_hat = I1I0_to_kappa(I1I0_th2_hat);
                    I1I0_th3_hat = sqrt( sum(Wc.*cos(SigPts(3,:)),2)^2 + sum(Wc.*sin(SigPts(3,:)),2)^2 );
                    kappa_th3_hat = I1I0_to_kappa(I1I0_th3_hat);
                case 'ut_standard'
                    P = zeros(dim_state,dim_state);
                    for j = 1:2*dim_state_aug+1
                        P = P + Wc(j)*SigPts(1:dim_state,j)*SigPts(1:dim_state,j)';
                    end
                    I1I0_th1_hat = sqrt( sum(Wc.*cos(SigPts(1,:)),2)^2 + sum(Wc.*sin(SigPts(1,:)),2)^2 );
                    kappa_th1_hat = I1I0_to_kappa(I1I0_th1_hat);
                    I1I0_th2_hat = sqrt( sum(Wc.*cos(SigPts(2,:)),2)^2 + sum(Wc.*sin(SigPts(2,:)),2)^2 );
                    kappa_th2_hat = I1I0_to_kappa(I1I0_th2_hat);
                    I1I0_th3_hat = sqrt( sum(Wc.*cos(SigPts(3,:)),2)^2 + sum(Wc.*sin(SigPts(3,:)),2)^2 );
                    kappa_th3_hat = I1I0_to_kappa(I1I0_th3_hat); 
                otherwise
                    error('UT type undefined');
            end

            %-------------------------- CORRECTION ----------------------------
            while obs(k) == bool
                switch sensor
                    case 'GPS'
                        ym = trajGT.trans(1:dimy,k) + sqrtR_true*randn(rstream,dimy,1);
                        Rk(1:dimy,1:dimy,k) = R;
                        Hk(1:dimy,1:dimy,k) = eye(dimy,dimy);
                        dim_obs_k(k) = dimy;
                    case 'Lidar'
                        [ym,index_visible_lm,dimy,sqrtR,sqrtR_true_tmp] = simul_meas_hlidar(trajGT.chi(:,:,k),landmarks,R_true,R,rstream);
                        if isempty(ym)
                            break
                        end
                        Rk(1:dimy,1:dimy,k) = sqrtR_true_tmp'*sqrtR_true_tmp;
                        H = (-1) * grad_h_lidar(Rot_GT,x_GT,landmarks,index_visible_lm);
                        Hk(1:dimy,1:6,k) = H;
                        dim_obs_k(k) = dimy;
                    otherwise
                        error('no sensor')
                end

                % re-generate sigma-points
                switch ut_type
                    case 'ut_sr'
                        S_aug = blkdiag(S,sqrtR); dim_state_aug = size(S_aug,1);
                        x_aug = zeros(dim_state_aug,1);
                        UT_kappa = 3-dim_state_aug;
                        [Wm,Wc,ksi,KSI,lambda] = compute_weights(dim_state_aug,UT_alpha,UT_beta,UT_kappa);
                        SigPts_01 = [zeros(dim_state_aug,1) -KSI*eye(dim_state_aug,dim_state_aug) KSI*eye(dim_state_aug,dim_state_aug)];
                        SigPts = zeros(dim_state_aug,2*dim_state_aug+1);
                        for i = 1:2*dim_state_aug+1
                            SigPts(:,i) = x_aug(:) + S_aug'*SigPts_01(:,i); % transpose if lower cholesky
                        end
                    case 'ut_gvm'
                        S_aug = blkdiag(S,sqrtR); dim_state_aug = size(S_aug,1);
                        state_aug = zeros(dim_state_aug,1);
                        n_aug = dim_state_aug-m; 
                        m_aug = dim_state_aug-n_aug; 
                        UT_kappa = 3-n_aug;
                        concentration_th1 = var_to_kappa(P(1,1)); concentration_th2 = var_to_kappa(P(2,2)); concentration_th3 = var_to_kappa(P(3,3));
                        switch correl_correction
                            case 'y'
                                Stmp = S;
                                M = P(4:dim_state,1:3)*P(1:3,1:3)^(-1);
                                Sxx_tilde = S_aug(4:dim_state,4:dim_state);
                                U = M*Stmp(1:3,1:3)'; % => U*U' = M*Pth*M';
                                %-------- replace Pxx = Pxx + MPthM' ----------
                                for j = 1:3
                                    Sxx_tilde = cholupdate(Sxx_tilde,U(:,j),'+'); % --->>>> Equivalent to Sxx_tilde = chol(P(4:dim_state,4:dim_state) + M*P(1:3,1:3)*M','upper') but more stable
                                end
                                %----------------------------------------------
                                S_aug(4:dim_state,4:dim_state) = Sxx_tilde;
                            case 'n'
                                %M = P(4:dim_state,1:3)*P(1:3,1:3)^(-1);
                                %S_aug(4:dim_state,4:dim_state) = chol(P(4:dim_state,4:dim_state)+M*P(1:3,1:3)*M','upper');
                        end
                        [Wm,Wc,ksi_x,ksi_th1,ksi_th2,ksi_th3,KSI] = compute_weights_vMG_update(n_aug,m_aug,UT_alpha,UT_kappa,concentration_th1,concentration_th2,concentration_th3); % Horwood Weights
                        SigPts_01 = [zeros(dim_state_aug,1) -KSI*eye(dim_state_aug,dim_state_aug) KSI*eye(dim_state_aug,dim_state_aug)];
                        SigPts = zeros(dim_state_aug,2*dim_state_aug+1);
                        for i = 1:2*dim_state_aug+1
                            SigPts(1:m,i)     = SigPts_01(1:m,i);
                            SigPts(m+1:end,i) = S_aug(m+1:end,m+1:end)'*SigPts_01(m+1:end,i); % do the transpose if "lower" cholesky
                        end
                    case 'ut_standard'
                        P_aug = blkdiag(P,sqrtR'*sqrtR); dim_state_aug = size(P_aug,1);
                        state_aug = zeros(dim_state_aug,1);
                        UT_kappa = 3-dim_state_aug;
                        [Wm,Wc,ksi,KSI,lambda] = compute_weights(dim_state_aug,UT_alpha,UT_beta,UT_kappa);
                        SigPts_01 = [zeros(dim_state_aug,1) -KSI*eye(dim_state_aug,dim_state_aug) KSI*eye(dim_state_aug,dim_state_aug)];
                        SigPts = zeros(dim_state_aug,2*dim_state_aug+1);
                        for i = 1:2*dim_state_aug+1
                            SigPts(:,i) = state_aug(:) + sqrtm(P_aug)*SigPts_01(:,i);
                        end
                    otherwise
                        error('UT type undefined');
                end

                switch CLVQbis
                    case 'yes'
                        gamma = diag([gamma_th1_bis;gamma_th2_bis;gamma_th3_bis;gamma_x_bis;gamma_y_bis;gamma_z_bis])*S; P = S'*S;
                        [SigPts,w_QO] = QO(gamma(ndimbis,ndimbis),Nbis,x_aug,P,SigPts,ndimbis);
                    case 'no'
                        ;
                    otherwise
                        error("CLVQ -> you have to specify 'yes' or 'no' ")
                end

                Y = zeros(dimy,2*dim_state_aug+1);
                for j = 1:2*dim_state_aug+1
                    ksi_j = SigPts(1:dim_state,j);
                    rj = SigPts(dim_state+1:end,j);
                    switch filter_side
                        case 'LEFT'
                            chi_j = chi*exp_multiSE3(ksi_j);
                        case 'RIGHT'
                            chi_j = exp_multiSE3(ksi_j)*chi;
                        otherwise
                            error('filter does not exist')
                    end
                    switch sensor
                        case 'GPS'
                            Y(:,j) = h_GPS(chi_j,rj);
                        case 'Lidar'
                            Y(:,j) = h_lidar(chi_j,landmarks,index_visible_lm,rj);
                        otherwise
                            error('no sensor')
                    end
                end

                % measurement prediction ypred, Py 
                switch ut_type
                    case 'ut_sr'
                        ypred = sum(Wm.*Y,2);
                        WY = sqrt(Wc(2:end)).*(Y(:,2:end)-ypred);
                        [~,RSy] = qr(WY');
                        Sy = RSy(1:dimy,1:dimy);
                        Uy = sqrt(abs(Wc(1)))*(Y(:,1) - ypred);
                        [Sy,~] = cholupdate(Sy,Uy,'-'); Py = Sy'*Sy;
                    case 'ut_gvm'
                        ypred = sum(Wm.*Y,2);
                        WY = sqrt(Wc(2:end)).*(Y(:,2:end)-ypred);
                        [~,RSy] = qr(WY');
                        Sy = RSy(1:dimy,1:dimy);
                        Uy = sqrt(abs(Wc(1)))*(Y(:,1) - ypred);
                        [Sy,~] = cholupdate(Sy,Uy,'-'); Py = Sy'*Sy;
                    case 'ut_standard'
                        ypred = sum(Wm.*Y,2);
                        Py = zeros(dimy,dimy);
                        for j = 1:2*dim_state_aug+1
                            Py = Py + Wc(j)*(Y(:,j)-ypred)*(Y(:,j)-ypred)';
                        end
                    otherwise
                        error('UT type undefined');
                end

                % cross-covariance Pxy
                Pxy = zeros(dim_state,dimy);
                for j = 2:2*dim_state_aug+1
                    Pxy = Pxy + Wc(j)*(SigPts(1:dim_state,j))*(Y(:,j)-ypred)';
                end

                % Kalman gain
                K = Pxy/Py;

                %-------------------------- update ----------------------------
                % innovation
                ksi_bar = K*(ym-ypred);

                % state update
                switch filter_side
                    case 'LEFT'
                        chi = chi*exp_multiSE3(ksi_bar);
                    case 'RIGHT'
                        chi = exp_multiSE3(ksi_bar)*chi;
                    otherwise
                        error('filter does not exist')
                end

                % covariance update
                switch ut_type
                    case 'ut_sr'
                        U = K*Sy';
                        for j = 1:dimy
                            [S,~] = cholupdate(S,U(:,j),'-');
                        end
                        J = JacSE3(wedge_se3(ksi_bar),filter_side);
                        S = S*J'; % set S = S*J' if P=S'S (upper cholesky) // set S = J*S if P=SS' (lower cholesky)
                        P = S'*S;
                    case 'ut_gvm'
                        U = K*Sy';
                        for j = 1:dimy
                            [S,~] = cholupdate(S,U(:,j),'-');
                        end
                        J = JacSE3(wedge_se3(ksi_bar),filter_side);
                        S = S*J'; % set S = S*J' if P=S'S (upper cholesky) // set S = J*S if P=SS' (lower cholesky)
                        P = S'*S;
                    case 'ut_standard'
                        P = P - K*Py*K';
                        J = JacSE3(wedge_se3(ksi_bar),filter_side);
                        P = J*P*J';
                    otherwise
                        error('UT type undefined');
                end

                break
            end % end correction

            % save data-----------------------------------------------
            Covs(:,:,k) = P;
            I1I0_hat_log(:,k) = [I1I0_th1_hat;I1I0_th2_hat;I1I0_th3_hat];
            kappa_hat_log(:,k) = [kappa_th1_hat;kappa_th2_hat;kappa_th3_hat];
            [Rot,eul_angles,x,add_pos] = chi2stateSE3(chi);
            trajukfLG = updateTrajSE3(trajukfLG,chi,Rot,eul_angles,x,k);
            LG_Error_k(:,k) = logSE3(invSE3(chiGT)*chi);
            %-----------------------------------------------------------
            
            T(k) = toc(start_filter);
        end %_________________________end filtering________________________
        delay = sum(T);

        % save data-----------------------------------------
        MC_computation_time(mc) = delay;
        MC_GT(mc).eul_angles = trajGT.eul_angles; MC_GT(mc).trans = trajGT.trans; MC_GT(mc).Rot = trajGT.Rot; MC_GT(mc).chi = trajGT.chi;
        MC_results(mc).eul_angles = trajukfLG.eul_angles; MC_results(mc).trans = trajukfLG.trans; MC_results(mc).Rot = trajukfLG.Rot; MC_results(mc).chi = trajukfLG.chi;
        MC_Covs(:,:,:,mc) = Covs;
        MC_I1I0_hat(:,:,mc) = I1I0_hat_log;
        MC_kappa_hat(:,:,mc) = kappa_hat_log;
        error_ukfLG = computeErrorSE3(trajukfLG,trajGT,k);
        MC_rmse_rot(mc) = sqrt(mean(error_ukfLG.errorRot.^2));
        MC_rmse_pos(mc) = sqrt(mean(error_ukfLG.errorPos.^2));
        MC_LG_Error_squared_k(:,:,mc) = LG_Error_k(:,:).^2;
        %-----------------------------------------------------

    end %_______________end monte-carlo simulation__________________________

    % save data for each noise level
    Noise_Level_MC_results(level,1:Ns) = MC_results; 
    Noise_Level_MC_GT(level,1:Ns) = MC_GT;
    Noise_Level_MC_Covs(:,:,:,:,level) = MC_Covs;
    Noise_Level_MC_I1I0_hat(:,:,:,level) = MC_I1I0_hat;

end %________________end sigma obs levels________________________________
kfin = k;

save('tmp.mat')

disp(['mean of the rotation RMSE over ' num2str(mc) ' monte-carlo simulations = ' num2str(mean(MC_rmse_rot)) '°'])
disp(['mean of the position RMSE over ' num2str(mc) ' monte-carlo simulations = ' num2str(mean(MC_rmse_pos)) 'm'])

figure(222);
plot3(trajGT.trans(1,1:kfin),trajGT.trans(2,1:kfin),trajGT.trans(3,1:kfin),'k','LineWidth',1); grid on; hold on;
plot3(trajukfLG.trans(1,1:kfin),trajukfLG.trans(2,1:kfin),trajukfLG.trans(3,1:kfin),'r','LineWidth',1);
plot3(trajGT.trans(1,1:1),trajGT.trans(2,1:1),trajGT.trans(3,1:1),'ok','LineWidth',1,'MarkerSize',10);
xlabel('$x_1$ (m)','interpreter','latex','FontSize',15)
ylabel('$x_2$ (m)','interpreter','latex','FontSize',15)
zlabel('$x_3$ (m)','interpreter','latex','FontSize',15) 
legend('Ground truth','Estimate','interpreter','latex','FontSize',14);

figure(111);
subplot(121)
plot(time(1:kfin),error_ukfLG.errorRot(1:kfin),'LineWidth',1); hold on; grid on;
xlabel('Time (sec)','interpreter','latex','FontSize',14); 
ylabel('Norm of the log-rotation error','interpreter','latex','FontSize',14)
subplot(122)
plot(time(1:kfin),error_ukfLG.errorPos(1:kfin),'LineWidth',1); hold on; grid on;
xlabel('Time (sec)','interpreter','latex','FontSize',14); 
ylabel('Norm of the log-position error','interpreter','latex','FontSize',14)
