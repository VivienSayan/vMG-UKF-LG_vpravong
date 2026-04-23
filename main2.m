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
f_obs = 5 % measurement frequency (Hz)

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
for level = 1:1
    level

    % Allocate data logs---------------
    MC_GT = struct('eul_angles',cell(1,Ns),'trans',cell(1,Ns),'Rot',cell(1,Ns),'chi',cell(1,Ns));
    
    MC_results_ukfLG = struct('eul_angles',cell(1,Ns),'trans',cell(1,Ns),'Rot',cell(1,Ns),'chi',cell(1,Ns));
    MC_Covs_ukfLG = zeros(6,6,kmax,Ns);
    MC_rmse_rot_ukfLG = zeros(1,Ns);
    MC_rmse_pos_ukfLG = zeros(1,Ns);
    MC_I1I0_hat_ukfLG = zeros(3,kmax,Ns);
    MC_kappa_hat_ukfLG = zeros(3,kmax,Ns);
    MC_LG_Error_squared_k_ukfLG = zeros(6,kmax,Ns);

    MC_results_vmgukfLG = struct('eul_angles',cell(1,Ns),'trans',cell(1,Ns),'Rot',cell(1,Ns),'chi',cell(1,Ns));
    MC_Covs_vmgukfLG = zeros(6,6,kmax,Ns);
    MC_rmse_rot_vmgukfLG = zeros(1,Ns);
    MC_rmse_pos_vmgukfLG = zeros(1,Ns);
    MC_I1I0_hat_vmgukfLG = zeros(3,kmax,Ns);
    MC_kappa_hat_vmgukfLG = zeros(3,kmax,Ns);
    MC_LG_Error_squared_k_vmgukfLG = zeros(6,kmax,Ns);
    %---------------------------------

    for mc = 1:Ns
        % initial state (true)
        eul_angles_GT = trajGT.eul_angles(:,1);
        x_GT = trajGT.trans(:,1);
        chiGT = state2chiSE3(eul_angles_GT,x_GT);
        [Rot_GT,eul_angles_GT,x_GT,add_pos] = chi2stateSE3(chiGT);
        trajGT = updateTrajSE3(trajGT,chiGT,Rot_GT,eul_angles_GT,x_GT,1);

        % initial covariance of the estimation
        p0Roll_ukfLG = (60*pi/180)^2; p0Pitch_ukfLG = (60*pi/180)^2; p0Yaw_ukfLG = (60*pi/180)^2; p0Rot_ukfLG = [p0Roll_ukfLG;p0Pitch_ukfLG;p0Yaw_ukfLG];
        p0Roll_vmgukfLG = p0Roll_ukfLG; p0Pitch_vmgukfLG = p0Pitch_ukfLG; p0Yaw_vmgukfLG = p0Yaw_ukfLG; p0Rot_vmgukfLG = [p0Roll_vmgukfLG;p0Pitch_vmgukfLG;p0Yaw_vmgukfLG];
        kappa_th1_hat_ukfLG = var_to_kappa(p0Roll_ukfLG); kappa_th2_hat_ukfLG = var_to_kappa(p0Pitch_ukfLG); kappa_th3_hat_ukfLG = var_to_kappa(p0Yaw_ukfLG);
        kappa_th1_hat_vmgukfLG = var_to_kappa(p0Roll_vmgukfLG); kappa_th2_hat_vmgukfLG = var_to_kappa(p0Pitch_vmgukfLG); kappa_th3_hat_vmgukfLG = var_to_kappa(p0Yaw_vmgukfLG);
        I1I0_th1_hat_ukfLG = bessel_ratio(0,kappa_th1_hat_ukfLG,10); I1I0_th2_hat_ukfLG = bessel_ratio(0,kappa_th2_hat_ukfLG,10); I1I0_th3_hat_ukfLG = bessel_ratio(0,kappa_th3_hat_ukfLG,10);
        I1I0_th1_hat_vmgukfLG = bessel_ratio(0,kappa_th1_hat_vmgukfLG,10); I1I0_th2_hat_vmgukfLG = bessel_ratio(0,kappa_th2_hat_vmgukfLG,10); I1I0_th3_hat_vmgukfLG = bessel_ratio(0,kappa_th3_hat_vmgukfLG,10);
        p0x_ukfLG = (0.5)^2; p0y_ukfLG = (0.5)^2; p0z_ukfLG = (0.5)^2; p0Pos_ukfLG = [p0x_ukfLG;p0y_ukfLG;p0z_ukfLG];
        p0x_vmgukfLG = (0.5)^2; p0y_vmgukfLG = (0.5)^2; p0z_vmgukfLG = (0.5)^2; p0Pos_vmgukfLG = [p0x_vmgukfLG;p0y_vmgukfLG;p0z_vmgukfLG];
        P_ukfLG = diag([p0Rot_ukfLG;p0Pos_ukfLG]); S_ukfLG = chol(P_ukfLG,'upper'); dim_state = size(P_ukfLG,1);
        P_vmgukfLG = diag([p0Rot_vmgukfLG;p0Pos_vmgukfLG]); S_vmgukfLG = chol(P_vmgukfLG,'upper');

        % initial estimate
        roll0_ukfLG = eul_angles_GT(3) + circ_vmrnd(0,kappa_th1_hat_ukfLG,1);%sqrt(p0Roll_ukfLG)*randn(rstream,1,1);% 
        roll0_vmgukfLG = roll0_ukfLG;
        pitch0_ukfLG = eul_angles_GT(2) + circ_vmrnd(0,kappa_th2_hat_ukfLG,1);%sqrt(p0Pitch_ukfLG)*randn(rstream,1,1);% 
        pitch0_vmgukfLG = pitch0_ukfLG;
        yaw0_ukfLG = eul_angles_GT(1) + circ_vmrnd(0,kappa_th3_hat_ukfLG,1);%sqrt(p0Yaw_ukfLG)*randn(rstream,1,1);% 
        yaw0_vmgukfLG = yaw0_ukfLG;
        eul_angles_ukfLG = [yaw0_ukfLG;pitch0_ukfLG;roll0_ukfLG];
        eul_angles_vmgukfLG = eul_angles_ukfLG;
        x_ukfLG = x_GT + sqrtm(diag(p0Pos_ukfLG))*randn(rstream,size(p0Pos_ukfLG,1),1);
        x_vmgukfLG = x_ukfLG;
        chi_ukfLG = state2chiSE3(eul_angles_ukfLG,x_ukfLG);
        chi_vmgukfLG = chi_ukfLG;

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
                r_x_true = sigma_model(1)*(0.01)^2; r_y_true = r_x_true; r_z_true = r_x_true; R_true = diag([r_x_true;r_y_true;r_z_true]);
                r_x = sigma_model(1)*(0.01)^2; r_y = r_x; r_z = r_x; R = diag([r_x;r_y;r_z]); 
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
        traj_ukfLG = initialTrajSE3(eul_angles_ukfLG,x_ukfLG,kmax);
        traj_vmgukfLG = initialTrajSE3(eul_angles_vmgukfLG,x_vmgukfLG,kmax);
        Covs_ukfLG = zeros(dim_state,dim_state,kmax); Covs_ukfLG(:,:,1) = S_ukfLG'*S_ukfLG;
        I1I0_hat_log_ukfLG = zeros(3,kmax); I1I0_hat_log_ukfLG(:,1) = [I1I0_th1_hat_ukfLG;I1I0_th2_hat_ukfLG;I1I0_th3_hat_ukfLG];
        kappa_hat_log_ukfLG = zeros(3,kmax); kappa_hat_log_ukfLG(:,1) = [kappa_th1_hat_ukfLG;kappa_th2_hat_ukfLG;kappa_th3_hat_ukfLG];
        MC_Covs_ukfLG(:,:,1,mc) = S_ukfLG'*S_ukfLG;
        MC_I1I0_hat_ukfLG(:,1,mc) = [I1I0_th1_hat_ukfLG;I1I0_th2_hat_ukfLG;I1I0_th3_hat_ukfLG];
        MC_kappa_hat_ukfLG(:,1,mc) = [kappa_th1_hat_ukfLG;kappa_th2_hat_ukfLG;kappa_th3_hat_ukfLG];
        LG_Error_k_ukfLG = zeros(6,kmax); LG_Error_k_ukfLG(:,1) = logSE3(invSE3(chiGT)*chi_ukfLG);

        Covs_vmgukfLG = zeros(dim_state,dim_state,kmax); Covs_vmgukfLG(:,:,1) = S_vmgukfLG'*S_vmgukfLG;
        I1I0_hat_log_vmgukfLG = zeros(3,kmax); I1I0_hat_log_vmgukfLG(:,1) = [I1I0_th1_hat_vmgukfLG;I1I0_th2_hat_vmgukfLG;I1I0_th3_hat_vmgukfLG];
        kappa_hat_log_vmgukfLG = zeros(3,kmax); kappa_hat_log_vmgukfLG(:,1) = [kappa_th1_hat_vmgukfLG;kappa_th2_hat_vmgukfLG;kappa_th3_hat_vmgukfLG];
        MC_Covs_vmgukfLG(:,:,1,mc) = S_vmgukfLG'*S_vmgukfLG;
        MC_I1I0_hat_vmgukfLG(:,1,mc) = [I1I0_th1_hat_vmgukfLG;I1I0_th2_hat_vmgukfLG;I1I0_th3_hat_vmgukfLG];
        MC_kappa_hat_vmgukfLG(:,1,mc) = [kappa_th1_hat_vmgukfLG;kappa_th2_hat_vmgukfLG;kappa_th3_hat_vmgukfLG];
        LG_Error_k_vmgukfLG = zeros(6,kmax); LG_Error_k_vmgukfLG(:,1) = logSE3(invSE3(chiGT)*chi_vmgukfLG);
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
            S_aug_ukfLG = blkdiag(S_ukfLG,sqrtQ); dim_state_aug = size(S_aug_ukfLG,1); % augment the state with the noises
            UT_kappa = 3-dim_state_aug;
            [Wm_ukfLG,Wc_ukfLG,ksi_ukfLG,KSI_ukfLG,lambda_ukfLG] = compute_weights(dim_state_aug,UT_alpha,UT_beta,UT_kappa); % Weights
            SigPts_01 = [zeros(dim_state_aug,1) -KSI_ukfLG*eye(dim_state_aug,dim_state_aug) KSI_ukfLG*eye(dim_state_aug,dim_state_aug)];
            SigPts_ukfLG = zeros(dim_state_aug,2*dim_state_aug+1);
            for i = 1:2*dim_state_aug+1
                SigPts_ukfLG(:,i) = S_aug_ukfLG'*SigPts_01(:,i); % do the transpose if "lower" cholesky
            end
            S_aug_vmgukfLG = blkdiag(S_vmgukfLG,sqrtQ); dim_state_aug = size(S_aug_vmgukfLG,1); % augment the state with the noises
            n_aug = dim_state_aug-2*m; % augmented gaussian variables
            m_aug = dim_state_aug-n_aug; % augmented circular variables
            UT_kappa = 3-n_aug;
            concentration_th1 = var_to_kappa(P_vmgukfLG(1,1)); concentration_th2 = var_to_kappa(P_vmgukfLG(2,2)); concentration_th3 = var_to_kappa(P_vmgukfLG(3,3));
            switch correl_prediction
                case 'y'
                    Stmp = S_vmgukfLG;
                    M = P_vmgukfLG(4:dim_state,1:3)*P_vmgukfLG(1:3,1:3)^(-1);
                    Sxx_tilde = S_aug_vmgukfLG(4:dim_state,4:dim_state);
                    U = M*Stmp(1:3,1:3)'; % => U*U' = M*Pth*M';
                    %-------- replace Pxx = Pxx + MPthM' ----------
                    for j = 1:3
                        Sxx_tilde = cholupdate(Sxx_tilde,U(:,j),'+'); % --->>>> Equivalent to Sxx_tilde = chol(P(4:dim_state,4:dim_state) + M*P(1:3,1:3)*M','upper') but more stable
                    end
                    %----------------------------------------------
                    S_aug_vmgukfLG(4:dim_state,4:dim_state) = Sxx_tilde;
                case 'n'
                    %M = P(4:dim_state,1:3)*P(1:3,1:3)^(-1);
                    %S_aug(4:dim_state,4:dim_state) = chol(P(4:dim_state,4:dim_state)+M*P(1:3,1:3)*M','upper');
            end
            [Wm_vmgukfLG,Wc_vmgukfLG,ksi_x,ksi_th1,ksi_th2,ksi_th3,KSI_vmgukfLG] = compute_weights_vMG_propagation(n_aug,m_aug,UT_alpha,UT_kappa,concentration_th1,concentration_th2,concentration_th3,kappa_q_th1,kappa_q_th2,kappa_q_th3); % Horwood Weights
            SigPts_01 = [zeros(dim_state_aug,1) -KSI_vmgukfLG*eye(dim_state_aug,dim_state_aug) KSI_vmgukfLG*eye(dim_state_aug,dim_state_aug)];
            SigPts_vmgukfLG = zeros(dim_state_aug,2*dim_state_aug+1);
            for i = 1:2*dim_state_aug+1
                SigPts_vmgukfLG(1:m,i)                     = SigPts_01(1:m,i);
                SigPts_vmgukfLG(m+1:dim_state,i)           = S_aug_vmgukfLG(m+1:dim_state,m+1:dim_state)'*SigPts_01(m+1:dim_state,i); % do the transpose if "lower" cholesky
                SigPts_vmgukfLG(dim_state+1:dim_state+m,i) = SigPts_01(dim_state+1:dim_state+m,i);
                SigPts_vmgukfLG(dim_state+m+1:end,i)       = S_aug_vmgukfLG(dim_state+m+1:end,dim_state+m+1:end)'*SigPts_01(dim_state+m+1:end,i);
            end

            switch CLVQ
                case 'yes'
                    gamma = diag([gamma_th1;gamma_th2;gamma_th3;gamma_x;gamma_y;gamma_z])*S_ukfLG; P_ukfLG = S_ukfLG'*S_ukfLG;
                    [SigPts_ukfLG,w_QO] = QO(gamma(ndim,ndim),N,x_aug,P_ukfLG,SigPts_ukfLG,ndim);
                case 'no'
                    ;
                otherwise
                    error("CLVQ -> you have to specify 'yes' or 'no' ")
            end

            chi_previous_ukfLG = chi_ukfLG; % save anterior pose
            chi_previous_vmgukfLG = chi_vmgukfLG;
            %-------------------------- PREDICTION ----------------------------        
            % mean propagation without noise
            vel = [u(1,k);u(2,k);u(3,k);u(4,k);u(5,k);u(6,k)];
            chi_ukfLG = chi_ukfLG*exp_multiSE3(vel*dt); % -> propagated mean state in the Lie group
            chi_vmgukfLG = chi_vmgukfLG*exp_multiSE3(vel*dt);
            chi_inv_ukfLG = invSE3(chi_ukfLG);
            chi_inv_vmgukfLG = invSE3(chi_vmgukfLG);

            for j = 2:2*dim_state_aug+1
                ksi_j_ukfLG = SigPts_ukfLG(1:dim_state,j);
                ksi_j_vmgukfLG = SigPts_vmgukfLG(1:dim_state,j);
                qj_ukfLG = SigPts_ukfLG(dim_state+1:end,j);
                qj_vmgukfLG = SigPts_vmgukfLG(dim_state+1:end,j);
                uj = vel;
                switch filter_side
                    case 'LEFT'
                        chi_j_ukfLG = chi_previous_ukfLG*exp_multiSE3(ksi_j_ukfLG);
                        chi_j_vmgukfLG = chi_previous_vmgukfLG*exp_multiSE3(ksi_j_vmgukfLG);
                        chi_j_ukfLG = chi_j_ukfLG*exp_multiSE3((uj+qj_ukfLG)*dt);
                        chi_j_vmgukfLG = chi_j_vmgukfLG*exp_multiSE3((uj+qj_vmgukfLG)*dt);
                        Xi_j_ukfLG = chi_inv_ukfLG*chi_j_ukfLG;
                        Xi_j_vmgukfLG = chi_inv_vmgukfLG*chi_j_vmgukfLG;
                    case 'RIGHT'
                        chi_j_ukfLG = exp_multiSE3(ksi_j_ukfLG)*chi_previous_ukfLG;
                        chi_j_vmgukfLG = exp_multiSE3(ksi_j_vmgukfLG)*chi_previous_vmgukfLG;
                        chi_j_ukfLG = chi_j_ukfLG*exp_multiSE3((uj+qj_ukfLG)*dt);
                        chi_j_vmgukfLG = chi_j_vmgukfLG*exp_multiSE3((uj+qj_vmgukfLG)*dt);
                        Xi_j_ukfLG = chi_j_ukfLG*chi_inv_ukfLG;
                        Xi_j_vmgukfLG = chi_j_vmgukfLG*chi_inv_vmgukfLG;
                    otherwise
                        error('filter does not exist');
                end
                SigPts_ukfLG(1:dim_state,j) = log_multiSE3(Xi_j_ukfLG); % propagated sigma-point j in the Lie algebra
                SigPts_vmgukfLG(1:dim_state,j) = log_multiSE3(Xi_j_vmgukfLG);
            end

            % compute a priori estimate and covariance
            WSigPts_ukfLG = sqrt(Wc_ukfLG(2:end)).*SigPts_ukfLG(1:dim_state,2:end);
            [~,RSx_ukfLG] = qr(WSigPts_ukfLG');
            S_ukfLG = RSx_ukfLG(1:dim_state,1:dim_state);
            P_ukfLG = S_ukfLG'*S_ukfLG;
            I1I0_th1_hat_ukfLG = sqrt( sum(Wc_ukfLG.*cos(SigPts_ukfLG(1,:)),2)^2 + sum(Wc_ukfLG.*sin(SigPts_ukfLG(1,:)),2)^2 );
            kappa_th1_hat_ukfLG = I1I0_to_kappa(I1I0_th1_hat_ukfLG);
            I1I0_th2_hat_ukfLG = sqrt( sum(Wc_ukfLG.*cos(SigPts_ukfLG(2,:)),2)^2 + sum(Wc_ukfLG.*sin(SigPts_ukfLG(2,:)),2)^2 );
            kappa_th2_hat_ukfLG = I1I0_to_kappa(I1I0_th2_hat_ukfLG);
            I1I0_th3_hat_ukfLG = sqrt( sum(Wc_ukfLG.*cos(SigPts_ukfLG(3,:)),2)^2 + sum(Wc_ukfLG.*sin(SigPts_ukfLG(3,:)),2)^2 );
            kappa_th3_hat_ukfLG = I1I0_to_kappa(I1I0_th3_hat_ukfLG);

            WSigPts_vmgukfLG = sqrt(Wc_vmgukfLG(2:end)).*SigPts_vmgukfLG(1:dim_state,2:end);
            [~,RSx_vmgukfLG] = qr(WSigPts_vmgukfLG');
            S_vmgukfLG = RSx_vmgukfLG(1:dim_state,1:dim_state);
            P_vmgukfLG = S_vmgukfLG'*S_vmgukfLG;
            I1I0_th1_hat_vmgukfLG = sqrt( sum(Wc_vmgukfLG.*cos(SigPts_vmgukfLG(1,:)),2)^2 + sum(Wc_vmgukfLG.*sin(SigPts_vmgukfLG(1,:)),2)^2 );
            kappa_th1_hat_vmgukfLG = I1I0_to_kappa(I1I0_th1_hat_vmgukfLG);
            I1I0_th2_hat_vmgukfLG = sqrt( sum(Wc_vmgukfLG.*cos(SigPts_vmgukfLG(2,:)),2)^2 + sum(Wc_vmgukfLG.*sin(SigPts_vmgukfLG(2,:)),2)^2 );
            kappa_th2_hat_vmgukfLG = I1I0_to_kappa(I1I0_th2_hat_vmgukfLG);
            I1I0_th3_hat_vmgukfLG = sqrt( sum(Wc_vmgukfLG.*cos(SigPts_vmgukfLG(3,:)),2)^2 + sum(Wc_vmgukfLG.*sin(SigPts_vmgukfLG(3,:)),2)^2 );
            kappa_th3_hat_vmgukfLG = I1I0_to_kappa(I1I0_th3_hat_vmgukfLG);

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
                S_aug_ukfLG = blkdiag(S_ukfLG,sqrtR); dim_state_aug = size(S_aug_ukfLG,1);
                x_aug = zeros(dim_state_aug,1);
                UT_kappa = 3-dim_state_aug;
                [Wm_ukfLG,Wc_ukfLG,ksi_ukfLG,KSI_ukfLG,lambda_ukfLG] = compute_weights(dim_state_aug,UT_alpha,UT_beta,UT_kappa);
                SigPts_01 = [zeros(dim_state_aug,1) -KSI_ukfLG*eye(dim_state_aug,dim_state_aug) KSI_ukfLG*eye(dim_state_aug,dim_state_aug)];
                SigPts_ukfLG = zeros(dim_state_aug,2*dim_state_aug+1);
                for i = 1:2*dim_state_aug+1
                    SigPts_ukfLG(:,i) = x_aug(:) + S_aug_ukfLG'*SigPts_01(:,i); % transpose if lower cholesky
                end

                S_aug_vmgukfLG = blkdiag(S_vmgukfLG,sqrtR); dim_state_aug = size(S_aug_vmgukfLG,1);
                n_aug = dim_state_aug-m; 
                m_aug = dim_state_aug-n_aug; 
                UT_kappa = 3-n_aug;
                concentration_th1 = var_to_kappa(P_vmgukfLG(1,1)); concentration_th2 = var_to_kappa(P_vmgukfLG(2,2)); concentration_th3 = var_to_kappa(P_vmgukfLG(3,3));
                switch correl_correction
                    case 'y'
                        Stmp = S_vmgukfLG;
                        M = P_vmgukfLG(4:dim_state,1:3)*P_vmgukfLG(1:3,1:3)^(-1);
                        Sxx_tilde = S_aug_vmgukfLG(4:dim_state,4:dim_state);
                        U = M*Stmp(1:3,1:3)'; % => U*U' = M*Pth*M';
                        %-------- replace Pxx = Pxx + MPthM' ----------
                        for j = 1:3
                            Sxx_tilde = cholupdate(Sxx_tilde,U(:,j),'+'); % --->>>> Equivalent to Sxx_tilde = chol(P(4:dim_state,4:dim_state) + M*P(1:3,1:3)*M','upper') but more stable
                        end
                        %----------------------------------------------
                        S_aug_vmgukfLG(4:dim_state,4:dim_state) = Sxx_tilde;
                    case 'n'
                        %M = P(4:dim_state,1:3)*P(1:3,1:3)^(-1);
                        %S_aug(4:dim_state,4:dim_state) = chol(P(4:dim_state,4:dim_state)+M*P(1:3,1:3)*M','upper');
                end
                [Wm_vmgukfLG,Wc_vmgukfLG,ksi_x,ksi_th1,ksi_th2,ksi_th3,KSI_vmgukfLG] = compute_weights_vMG_update(n_aug,m_aug,UT_alpha,UT_kappa,concentration_th1,concentration_th2,concentration_th3); % Horwood Weights
                SigPts_01 = [zeros(dim_state_aug,1) -KSI_vmgukfLG*eye(dim_state_aug,dim_state_aug) KSI_vmgukfLG*eye(dim_state_aug,dim_state_aug)];
                SigPts_vmgukfLG = zeros(dim_state_aug,2*dim_state_aug+1);
                for i = 1:2*dim_state_aug+1
                    SigPts_vmgukfLG(1:m,i)     = SigPts_01(1:m,i);
                    SigPts_vmgukfLG(m+1:end,i) = S_aug_vmgukfLG(m+1:end,m+1:end)'*SigPts_01(m+1:end,i); % do the transpose if "lower" cholesky
                end

                switch CLVQbis
                    case 'yes'
                        gamma = diag([gamma_th1_bis;gamma_th2_bis;gamma_th3_bis;gamma_x_bis;gamma_y_bis;gamma_z_bis])*S_ukfLG; P_ukfLG = S_ukfLG'*S_ukfLG;
                        [SigPts_ukfLG,w_QO] = QO(gamma(ndimbis,ndimbis),Nbis,x_aug,P_ukfLG,SigPts_ukfLG,ndimbis);
                    case 'no'
                        ;
                    otherwise
                        error("CLVQ -> you have to specify 'yes' or 'no' ")
                end

                Y_ukfLG = zeros(dimy,2*dim_state_aug+1);
                Y_vmgukfLG = Y_ukfLG;
                for j = 1:2*dim_state_aug+1
                    ksi_j_ukfLG = SigPts_ukfLG(1:dim_state,j);
                    ksi_j_vmgukfLG = SigPts_vmgukfLG(1:dim_state,j);
                    rj_ukfLG = SigPts_ukfLG(dim_state+1:end,j);
                    rj_vmgukfLG = SigPts_vmgukfLG(dim_state+1:end,j);
                    switch filter_side
                        case 'LEFT'
                            chi_j_ukfLG = chi_ukfLG*exp_multiSE3(ksi_j_ukfLG);
                            chi_j_vmgukfLG = chi_vmgukfLG*exp_multiSE3(ksi_j_vmgukfLG);
                        case 'RIGHT'
                            chi_j_ukfLG = exp_multiSE3(ksi_j_ukfLG)*chi_ukfLG;
                            chi_j_vmgukfLG = exp_multiSE3(ksi_j_vmgukfLG)*chi_vmgukfLG;
                        otherwise
                            error('filter does not exist')
                    end
                    switch sensor
                        case 'GPS'
                            Y_ukfLG(:,j) = h_GPS(chi_j_ukfLG,rj_ukfLG);
                            Y_vmgukfLG(:,j) = h_GPS(chi_j_vmgukfLG,rj_vmgukfLG);
                        case 'Lidar'
                            Y_ukfLG(:,j) = h_lidar(chi_j_ukfLG,landmarks,index_visible_lm,rj_ukfLG);
                            Y_vmgukfLG(:,j) = h_lidar(chi_j_vmgukfLG,landmarks,index_visible_lm,rj_vmgukfLG);
                        otherwise
                            error('no sensor')
                    end
                end

                % measurement prediction ypred, Py 
                ypred_ukfLG = sum(Wm_ukfLG.*Y_ukfLG,2);
                WY_ukfLG = sqrt(Wc_ukfLG(2:end)).*(Y_ukfLG(:,2:end)-ypred_ukfLG);
                [~,RSy_ukfLG] = qr(WY_ukfLG');
                Sy_ukfLG = RSy_ukfLG(1:dimy,1:dimy);
                Uy = sqrt(abs(Wc_ukfLG(1)))*(Y_ukfLG(:,1) - ypred_ukfLG);
                [Sy_ukfLG,~] = cholupdate(Sy_ukfLG,Uy,'-'); Py_ukfLG = Sy_ukfLG'*Sy_ukfLG;

                ypred_vmgukfLG = sum(Wm_vmgukfLG.*Y_vmgukfLG,2);
                WY_vmgukfLG = sqrt(Wc_vmgukfLG(2:end)).*(Y_vmgukfLG(:,2:end)-ypred_vmgukfLG);
                [~,RSy_vmgukfLG] = qr(WY_vmgukfLG');
                Sy_vmgukfLG = RSy_vmgukfLG(1:dimy,1:dimy);
                Uy = sqrt(abs(Wc_vmgukfLG(1)))*(Y_vmgukfLG(:,1) - ypred_vmgukfLG);
                [Sy_vmgukfLG,~] = cholupdate(Sy_vmgukfLG,Uy,'-'); Py_vmgukfLG = Sy_vmgukfLG'*Sy_vmgukfLG;

                % cross-covariance Pxy
                Pxy_ukfLG = zeros(dim_state,dimy);
                Pxy_vmgukfLG = Pxy_ukfLG;
                for j = 2:2*dim_state_aug+1
                    Pxy_ukfLG = Pxy_ukfLG + Wc_ukfLG(j)*(SigPts_ukfLG(1:dim_state,j))*(Y_ukfLG(:,j)-ypred_ukfLG)';
                    Pxy_vmgukfLG = Pxy_vmgukfLG + Wc_vmgukfLG(j)*(SigPts_vmgukfLG(1:dim_state,j))*(Y_vmgukfLG(:,j)-ypred_vmgukfLG)';
                end

                % Kalman gain
                K_ukfLG = Pxy_ukfLG/Py_ukfLG;
                K_vmgukfLG = Pxy_vmgukfLG/Py_vmgukfLG;

                %-------------------------- update ----------------------------
                % innovation
                ksi_bar_ukfLG = K_ukfLG*(ym-ypred_ukfLG);
                ksi_bar_vmgukfLG = K_vmgukfLG*(ym-ypred_vmgukfLG);

                % state update
                switch filter_side
                    case 'LEFT'
                        chi_ukfLG = chi_ukfLG*exp_multiSE3(ksi_bar_ukfLG);
                        chi_vmgukfLG = chi_vmgukfLG*exp_multiSE3(ksi_bar_vmgukfLG);
                    case 'RIGHT'
                        chi_ukfLG = exp_multiSE3(ksi_bar_ukfLG)*chi_ukfLG;
                        chi_vmgukfLG = exp_multiSE3(ksi_bar_vmgukfLG)*chi_vmgukfLG;
                    otherwise
                        error('filter does not exist')
                end

                % covariance update
                U = K_ukfLG*Sy_ukfLG';
                for j = 1:dimy
                    [S_ukfLG,~] = cholupdate(S_ukfLG,U(:,j),'-');
                end
                J = JacSE3(wedge_se3(ksi_bar_ukfLG),filter_side);
                S_ukfLG = S_ukfLG*J'; % set S = S*J' if P=S'S (upper cholesky) // set S = J*S if P=SS' (lower cholesky)
                P_ukfLG = S_ukfLG'*S_ukfLG;

                U = K_vmgukfLG*Sy_vmgukfLG';
                for j = 1:dimy
                    [S_vmgukfLG,~] = cholupdate(S_vmgukfLG,U(:,j),'-');
                end
                J = JacSE3(wedge_se3(ksi_bar_vmgukfLG),filter_side);
                S_vmgukfLG = S_vmgukfLG*J'; % set S = S*J' if P=S'S (upper cholesky) // set S = J*S if P=SS' (lower cholesky)
                P_vmgukfLG = S_vmgukfLG'*S_vmgukfLG;

                break
            end % end correction

            % save data-----------------------------------------------
            [Rot_ukfLG,eul_angles_ukfLG,x_ukfLG,add_pos] = chi2stateSE3(chi_ukfLG);
            traj_ukfLG = updateTrajSE3(traj_ukfLG,chi_ukfLG,Rot_ukfLG,eul_angles_ukfLG,x_ukfLG,k);
            Covs_ukfLG(:,:,k) = P_ukfLG;
            I1I0_hat_log_ukfLG(:,k) = [I1I0_th1_hat_ukfLG;I1I0_th2_hat_ukfLG;I1I0_th3_hat_ukfLG];
            kappa_hat_log_ukfLG(:,k) = [kappa_th1_hat_ukfLG;kappa_th2_hat_ukfLG;kappa_th3_hat_ukfLG];
            LG_Error_k_ukfLG(:,k) = logSE3(invSE3(chiGT)*chi_ukfLG);

            [Rot_vmgukfLG,eul_angles_vmgukfLG,x_vmgukfLG,add_pos] = chi2stateSE3(chi_vmgukfLG);
            traj_vmgukfLG = updateTrajSE3(traj_vmgukfLG,chi_vmgukfLG,Rot_vmgukfLG,eul_angles_vmgukfLG,x_vmgukfLG,k);
            Covs_vmgukfLG(:,:,k) = P_vmgukfLG;
            I1I0_hat_log_vmgukfLG(:,k) = [I1I0_th1_hat_vmgukfLG;I1I0_th2_hat_vmgukfLG;I1I0_th3_hat_vmgukfLG];
            kappa_hat_log_vmgukfLG(:,k) = [kappa_th1_hat_vmgukfLG;kappa_th2_hat_vmgukfLG;kappa_th3_hat_vmgukfLG];
            LG_Error_k_vmgukfLG(:,k) = logSE3(invSE3(chiGT)*chi_vmgukfLG);
            %-----------------------------------------------------------
            
            T(k) = toc(start_filter);
        end %_________________________end filtering________________________
        delay = sum(T);

        % save data-----------------------------------------
        MC_GT(mc).eul_angles = trajGT.eul_angles; MC_GT(mc).trans = trajGT.trans; MC_GT(mc).Rot = trajGT.Rot; MC_GT(mc).chi = trajGT.chi;
        
        MC_results_ukfLG(mc).eul_angles = traj_ukfLG.eul_angles; MC_results_ukfLG(mc).trans = traj_ukfLG.trans; MC_results_ukfLG(mc).Rot = traj_ukfLG.Rot; MC_results_ukfLG(mc).chi = traj_ukfLG.chi;
        MC_Covs_ukfLG(:,:,:,mc) = Covs_ukfLG;
        MC_I1I0_hat_ukfLG(:,:,mc) = I1I0_hat_log_ukfLG;
        MC_kappa_hat_ukfLG(:,:,mc) = kappa_hat_log_ukfLG;
        error_ukfLG = computeErrorSE3(traj_ukfLG,trajGT,k);
        MC_rmse_rot_ukfLG(mc) = sqrt(mean(error_ukfLG.errorRot.^2));
        MC_rmse_pos_ukfLG(mc) = sqrt(mean(error_ukfLG.errorPos.^2));
        MC_LG_Error_squared_k_ukfLG(:,:,mc) = LG_Error_k_ukfLG(:,:).^2;

        MC_results_vmgukfLG(mc).eul_angles = traj_vmgukfLG.eul_angles; MC_results_vmgukfLG(mc).trans = traj_vmgukfLG.trans; MC_results_vmgukfLG(mc).Rot = traj_vmgukfLG.Rot; MC_results_vmgukfLG(mc).chi = traj_vmgukfLG.chi;
        MC_Covs_vmgukfLG(:,:,:,mc) = Covs_vmgukfLG;
        MC_I1I0_hat_vmgukfLG(:,:,mc) = I1I0_hat_log_vmgukfLG;
        MC_kappa_hat_vmgukfLG(:,:,mc) = kappa_hat_log_vmgukfLG;
        error_vmgukfLG = computeErrorSE3(traj_vmgukfLG,trajGT,k);
        MC_rmse_rot_vmgukfLG(mc) = sqrt(mean(error_vmgukfLG.errorRot.^2));
        MC_rmse_pos_vmgukfLG(mc) = sqrt(mean(error_vmgukfLG.errorPos.^2));
        MC_LG_Error_squared_k_vmgukfLG(:,:,mc) = LG_Error_k_vmgukfLG(:,:).^2;
        %-----------------------------------------------------

    end %_______________end monte-carlo simulation__________________________

    % save data for each noise level
    Noise_Level_MC_GT_ukfLG(level,1:Ns) = MC_GT;

    Noise_Level_MC_results_ukfLG(level,1:Ns) = MC_results_ukfLG; 
    Noise_Level_MC_Covs_ukfLG(:,:,:,:,level) = MC_Covs_ukfLG;
    Noise_Level_MC_I1I0_hat_ukfLG(:,:,:,level) = MC_I1I0_hat_ukfLG;

    Noise_Level_MC_results_vmgukfLG(level,1:Ns) = MC_results_vmgukfLG; 
    Noise_Level_MC_Covs_vmgukfLG(:,:,:,:,level) = MC_Covs_vmgukfLG;
    Noise_Level_MC_I1I0_hat_vmgukfLG(:,:,:,level) = MC_I1I0_hat_vmgukfLG;

end %________________end sigma obs levels________________________________
kfin = k;

save('tmp.mat')

disp(['UKF-LG: mean of the rotation RMSE over ' num2str(mc) ' monte-carlo simulations = ' num2str(mean(MC_rmse_rot_ukfLG)) '°'])
disp(['UKF-LG: mean of the position RMSE over ' num2str(mc) ' monte-carlo simulations = ' num2str(mean(MC_rmse_pos_ukfLG)) 'm'])
disp(' ');
disp(['vMG-UKF-LG: mean of the rotation RMSE over ' num2str(mc) ' monte-carlo simulations = ' num2str(mean(MC_rmse_rot_vmgukfLG)) '°'])
disp(['vMG-UKF-LG: mean of the position RMSE over ' num2str(mc) ' monte-carlo simulations = ' num2str(mean(MC_rmse_pos_vmgukfLG)) 'm'])

figure(222);
plot3(trajGT.trans(1,1:kfin),trajGT.trans(2,1:kfin),trajGT.trans(3,1:kfin),'k','LineWidth',1); grid on; hold on;
plot3(traj_ukfLG.trans(1,1:kfin),traj_ukfLG.trans(2,1:kfin),traj_ukfLG.trans(3,1:kfin),'LineWidth',1);
plot3(traj_vmgukfLG.trans(1,1:kfin),traj_vmgukfLG.trans(2,1:kfin),traj_vmgukfLG.trans(3,1:kfin),'LineWidth',1);
plot3(trajGT.trans(1,1:1),trajGT.trans(2,1:1),trajGT.trans(3,1:1),'ok','LineWidth',1,'MarkerSize',10);
xlabel('$x_1$ (m)','interpreter','latex','FontSize',15)
ylabel('$x_2$ (m)','interpreter','latex','FontSize',15)
zlabel('$x_3$ (m)','interpreter','latex','FontSize',15) 
legend('Ground truth','UKF-LG','vMG-UKF-LG','interpreter','latex','FontSize',14);

figure(111);
subplot(121)
plot(time(1:kfin),error_ukfLG.errorRot(1:kfin),'LineWidth',1); hold on; grid on;
plot(time(1:kfin),error_vmgukfLG.errorRot(1:kfin),'LineWidth',1);
xlabel('Time (sec)','interpreter','latex','FontSize',14); 
ylabel('Norm of the log-rotation error','interpreter','latex','FontSize',14)
legend('UKF-LG','vMG-UKF-LG','interpreter','latex','FontSize',14);
subplot(122)
plot(time(1:kfin),error_ukfLG.errorPos(1:kfin),'LineWidth',1); hold on; grid on;
plot(time(1:kfin),error_vmgukfLG.errorPos(1:kfin),'LineWidth',1);
xlabel('Time (sec)','interpreter','latex','FontSize',14); 
ylabel('Norm of the log-position error','interpreter','latex','FontSize',14)
legend('UKF-LG','vMG-UKF-LG','interpreter','latex','FontSize',14);
