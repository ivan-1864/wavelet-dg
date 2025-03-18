    % close all
clear

timetime = clock;

D=importdata('../../Data/Output_data/anomaly_data.txt');
true_anomaly_all=D.data;
clear D

D=importdata('../../Data/Output_data/anomaly_low_data.txt');
true_anomaly_low=D.data;
clear D

true_anomaly_res = [true_anomaly_low(:, 1), true_anomaly_all(:, 2:end) - true_anomaly_low(:, 2:end)];

D=importdata('../../Data/Output_data/anomaly_sat_data.txt');
true_anomaly_sat=D.data;
clear D

% --------------download data----------------------------------------------
D               = importdata('../../INS/Output_data/output_INS_data.txt');  
time            = D.data(:, 1);
dVx_array       = D.data(:, 2:4);
Vx_array        = D.data(:, 5:7);
Omega_x_array   = D.data(:, 8:10);
u_E_array       = D.data(:, 23:25);
L_zx_array      = [D.data(:, 11), D.data(:, 12), D.data(:, 13),...
                   D.data(:, 14), D.data(:, 15), D.data(:, 16),...
                   D.data(:, 17), D.data(:, 18), D.data(:, 19)];
g0_array        = D.data(:, 20:22); 
omega_x_array   = Omega_x_array + u_E_array;
clear D



cfg; % download configuration param.
% ---------------global variables--------------------------------------
Length          = size(time, 1);
time_step       = mean(diff(time));
% -------------------------------------------------------------------------


t_val           = time;
t_min           = min(t_val);
t_max           = max(t_val);
TimeWave = t_val - t_min; 

if TimeEnd == 0
    TimeEnd = Length;
end

% k_arr = [30, 60, 60];
k_arr = [30, 60];
j_min = -8;
% j_min = -9;
j_max = -7;

% j_max = -6;
len_psi = sum(k_arr)*3;

% ----------creating matrix P, Q, R (for KF)-------------------------------

d_dv            = Std_dV * ones(1, 3);
d_beta          = deg2rad(Std_beta / 60) * ones(1, 3);
d_dg            = Std_dg * 10^(-5) * ones(1, 3);
d_p             = Std_p * 10^(-5) * ones(1, 3);
d_f             = [Std_f * 10^(-5) * ones(1, 2), 10^(-5)];
d_nu            = deg2rad(Std_nu) / 3600 * ones(1, 3);
d_w = 1000*10^(-5)*ones(1, len_psi);
P_last_diag     = [d_dv, d_beta, d_f, d_nu, d_w]; % diagonal of sqrt of cov. matrix at time t=0
P_last          = diag(P_last_diag);
% -------------------------------------------------------------------------
d_ax            = (10^(-5) * Std_AX) * ones(1, 3);
d_dus           = (deg2rad(Std_DUS) / 3600) * ones(1, 3);
% d_q             = (10^(-5) * Std_q) * ones(1, 3);
Q_sqrt_diag     = [d_ax, d_dus];
Q_sqrt          = diag(Q_sqrt_diag);
% -------------------------------------------------------------------------
d_vel           = (Std_GPS) * ones(1, 3);
R_sqrt_diag     = d_vel;
R_sqrt          = diag(R_sqrt_diag);
R_sqrt_inv      = diag(1 ./ R_sqrt_diag);


% --------------------------state vector--------------------
dimX            = 12+len_psi;  %  x = (dlt V, dlt F, nu, beta, dg, p)
Y_last          = zeros(dimX, 1);
%----------------- Declare arrays for results ---------------------- 
Y_forw          = zeros(TimeEnd, dimX);
P_forw          = cell(1, TimeEnd);
Tmp_matr1       = cell(1,TimeEnd);               % Auxiliary matrix
Tmp_matr2       = cell(1,TimeEnd);               % Auxiliary matrix
Resid_forw      = zeros(TimeEnd,size(R_sqrt,1)); % Array of residuals

P_forw{1}       = P_last;

X_j=zeros(dimX, TimeEnd);

Psi = zeros(1, sum(k_arr));
% Psi3d = 
% ------------------KF algorithm-------------------------------------------
disp('начало работы алгоритма')

for i = 1:TimeEnd
%    --------------initialization-------------------------------------------
    for j = j_min : j_max
        for k = 1 : k_arr(j - j_min + 1)
           dT  = (2^j) * t_max / k_arr(j - j_min + 1);  % step of wavelet grid
           t_k = k * dT;                % knot in grid
           Psi(:, sum(k_arr(1: j - j_min)) + k) = wavel_trf(j,t_k,TimeWave(i));   
        end
    end

    Psi3d = [Psi, zeros(size(Psi)), zeros(size(Psi));
             zeros(size(Psi)),  Psi, zeros(size(Psi));
             zeros(size(Psi)), zeros(size(Psi)), Psi;];
    g           = Make_Cross_Matrix(g0_array(i, :)); % true_anomaly_sat(10*i, 2:end));
    Vx          = Make_Cross_Matrix(Vx_array(i, :));
    Z_t         = dVx_array(i, :)';
    omega_x     = Make_Cross_Matrix(omega_x_array(i, :));
    omega_pl_u  = Make_Cross_Matrix(Omega_x_array(i, :)) + 2 * Make_Cross_Matrix(u_E_array(i, :));
%     ---------------------------------------------------------------------
    L_zx        = [L_zx_array(i,1), L_zx_array(i,2), L_zx_array(i,3);
                   L_zx_array(i,4), L_zx_array(i,5), L_zx_array(i,6);
                   L_zx_array(i,7), L_zx_array(i,8), L_zx_array(i,9)];
    len_psi = length(Psi3d);
%     ---------------------------------------------------------------------
    A           = [omega_pl_u, g, L_zx' , Vx * L_zx', -Psi3d;
                   zeros(3), omega_x, zeros(3), L_zx', zeros(3, len_psi);
                   zeros(6, 12+len_psi);
                   zeros(len_psi, 12+len_psi)];
       
    H_t         = [eye(3), -Vx, zeros(3, 6), zeros(3, len_psi)];
    F_t         =  eye(size(A)) + A * time_step;
%     ---------------------------------------------------------------------
    J_t         = [L_zx', Vx * L_zx';
                  zeros(3), L_zx';         
                  zeros(6);
                  zeros(len_psi, 6)];
    
    [Y_pred,P_pred,Tmp1,Tmp2,Resid] = KF_forward( ...
        F_t,J_t,Q_sqrt,Z_t,H_t,R_sqrt,R_sqrt_inv,Y_last,P_last);
    Y_last = Y_pred; 
    P_last = P_pred;     
 
    X_j(:, i) = P_pred * Y_pred;

    if i/100 - floor(i/100) == 0
        disp(num2str(i));
    end
end



% ---------------------record into the file--------------------------------

time_end = 5250;

% % dv
% figure(1)
% hold
% title('dV')
% plot(time(1:time_end*10), 10^5*X_sm(1:end-1, 1:3))
% % 
% % beta
% figure(2)
% hold
% title('beta')
% plot(time(1:time_end*10), 10^5*X_sm(1:end-1, 4:6))
% % 
% % 
% % df
% figure(3)
% hold
% title('Delta f')
% plot(time(1:time_end*10), 10^5*X_sm(1:end-1, 7:9))
% 
% % nu
% figure(4)
% hold
% title('nu')
% plot(time(1:time_end*10), 10^5*X_sm(1:end-1, 10:12))

time_end = TimeEnd/10;


% dg3
% figure(5)
% hold
% plot(time(1:time_end*10), 10^5*X_j(15, :))
% plot(true_anomaly_res(1*10:10:time_end*100,1),10^5*true_anomaly_res(1:10:time_end*100,4));
% legend('comp', 'true')
% 
% figure(1)
% hold
% plot(time(1:time_end*10), 10^5*X_j(15, :))
% plot(true_anomaly_res(:, 1),10^5*true_anomaly_res(:,4), LineWidth=1.5);
% plot(true_anomaly_low(:, 1),10^5*true_anomaly_low(:,4), LineWidth=1.5);
% plot(true_anomaly_sat(:, 1),10^5*true_anomaly_all(:,4), LineWidth=1.5);
% legend('est', 'res', 'low', 'all')

% figure(6)
% hold
% plot(time(1:time_end*10), 10^5*X_j(13, :))
% plot(true_anomaly_res(1:10:time_end*100,1), 10^5*true_anomaly_res(1:10:time_end*100,2));
% legend('comp', 'true')
% 
% figure(7)
% hold
% plot(time(1:time_end*10), 10^5*X_j(14, :))
% plot(true_anomaly_res(1:10:time_end*100,1), 10^5*true_anomaly_res(1:10:time_end*100,3));
% legend('comp', 'true')

disp(['время работы программы: ', num2str(etime(clock, timetime)), ' секунд'])