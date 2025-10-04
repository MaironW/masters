% MP208
% Computational Exercise 5
% Mairon de Souza Wolniewicz
% Written in MATLAB R2018b

%% Setup
clc; clear; close all;

% Plot colors
c0 = [0.0,0.0,0.0];
c1 = [0.0,0.4,1.0];
c2 = [0.8,0.0,0.0];
c3 = [0.2,0.8,0.0];
c4 = [1.0,0.0,1.0];
c1w = c1 + 0.7*(1-c1);
c1c3 = (1-0.5)*c1 + 0.5*c3;
c2c4 = (1-0.5)*c4 + 0.5*c2;

% Simulation parameters
N_mc    =  100; % Number of realizations
t_start =  0.0; % [s] Start time
t_end   = 20.0; % [s] End time
t_step  =  0.1; % [s] Sample time (T from the problem statement)
t_range = t_start:t_step:t_end;
N_steps = length(t_range)-1; % Number of steps
N_states = 2;

% Initial values
x_bar = zeros(N_states,1); % Initial state mean
P_bar = eye(N_states);     % Initial state covariance

% Model parameters
G = eye(N_states);
Q = 0.01*eye(N_states); % Process covariance matrix
R = 0.01;               % Measurement variance

% Functions and Jacobians
f = @(x, u) [-x(1,:) + x(2,:); -0.1*x(1,:).^2 - 1 + u];
h = @(x) x(1,:);
F = @(x) [-1 1; -0.2*x(1) 0];
H = [1 0];
P_dot = @(x,F,P,G,Q) F(x)*P + P*F(x)' + G*Q*G';

% Preset plots so lines are generated in computation order
ax_state   = plot_setup({'System Output y(t)','Control Input u(t)','State x1(t)','State x2(t)'},{'','','','t (s)'},{'y(t)','u(t)','x1(t)','x2(t)'});
ax_error   = plot_setup({'State x1(t) Error','State x2(t) Error','Innovation'},{'','','t (s)'},{'x1(t)','x2(t)','Innov'});
ukf_vs_ekf = plot_setup({'State x1(t) Error','State x2(t) Error','Innovation'},{'','','t (s)'},{'x1(t)','x2(t)','Innov'});

%% Run Monter-Carlo realizations
y_mc      = zeros(N_steps+1, N_mc);
u_mc      = zeros(N_steps+1, N_mc);
x_mc      = zeros(2, N_steps+1, N_mc);
x_est_mc  = zeros(2, N_steps+1, N_mc);
innov_mc  = zeros(N_steps+1, N_mc);
P_diag_mc = zeros(2, N_steps+1, N_mc);

%% a) Simulate the system
for iter_mc = 1:N_mc
    x = zeros(2,N_steps+1);
    y = zeros(1,N_steps+1);
    u = zeros(1,N_steps+1);

    % Realizations
    w      = sqrt(Q) * randn(2,N_steps+1);
    v      = sqrt(R) * randn(1,N_steps+1);
    x(:,1) = x_bar + sqrt(P_bar) * randn(2,1);
    y(1)   = h(x(:,1)) + v(1);
    u(1)   = -10*y(1) + 10;

    % Integration
    for k = 1:N_steps
        x(:,k+1) = rk4_x(f,x(:,k),u(k),w(:,k),t_step);

        y(k+1) = h(x(:,k+1)) + v(k+1);
        u(k+1) = -10*y(k+1) + 10;
    end

    y_mc(:,iter_mc)   = y;
    u_mc(:,iter_mc)   = u;
    x_mc(:,:,iter_mc) = x;

    plot(ax_state(1),t_range,y,'Color',c1w,'LineWidth',1,'HandleVisibility','off');
    plot(ax_state(2),t_range,u,'Color',c1w,'LineWidth',1,'HandleVisibility','off');
    plot(ax_state(3),t_range,x(1,:)','Color',c1w,'LineWidth',1,'HandleVisibility','off');
    plot(ax_state(4),t_range,x(2,:)','Color',c1w,'LineWidth',1,'HandleVisibility','off');

end
% Add legends to Monte-Carlo iteration
mc_str = [num2str(N_mc),' MC runs'];
plot(ax_state(1),nan,nan,'Color',c1w,'LineWidth',1,'DisplayName',['y(t) ',mc_str]);
plot(ax_state(2),nan,nan,'Color',c1w,'LineWidth',1,'DisplayName',['u(t) ',mc_str]);
plot(ax_state(3),nan,nan,'Color',c1w,'LineWidth',1,'DisplayName',['x1(t) ',mc_str]);
plot(ax_state(4),nan,nan,'Color',c1w,'LineWidth',1,'DisplayName',['x2(t) ',mc_str]);

% Plot mean values of system signals
y_mean = mean(y_mc,2);
u_mean = mean(u_mc,2);
x_mean = mean(x_mc,3);
plot(ax_state(1),t_range,y_mean','Color',c1,'LineWidth',2,'DisplayName','y(t) Real Mean');
plot(ax_state(2),t_range,u_mean','Color',c1,'LineWidth',2,'DisplayName','u(t) Real Mean');
plot(ax_state(3),t_range,x_mean(1,:)','Color',c1,'LineWidth',2,'DisplayName','x1(t) Real Mean');
plot(ax_state(4),t_range,x_mean(2,:)','Color',c1,'LineWidth',2,'DisplayName','x2(t) Real Mean');

%% b) Design and implement a CDEnKF for estimating {x(t)}
tic;
for iter_mc = 1:N_mc
    [x_est, P, innov] = ensemble_kalman_filter(f,h,G,R,Q,u_mc(:,iter_mc),y_mc(:,iter_mc),t_step,x_bar,P_bar,100);
    x_est_mc(:,:,iter_mc) = x_est;
    innov_mc(:,iter_mc) = innov;
    P_diag_mc(:,:,iter_mc) = [P(1,1,:); P(2,2,:)];
end
EnKF_dt = toc;
disp(['EnKF time to run: ',num2str(EnKF_dt), ' s']);

%% c) Evaluate CDEnKF performance

% Calculate statistics
x_err_mc   = x_est_mc - x_mc;
x_err_mean = mean(x_err_mc,3);
x_est_mean = mean(x_est_mc,3);
innov_mean = mean(innov_mc,2);
x_est_rms  = sqrt(mean(x_err_mc.^2, 3));
P_std = sqrt(mean(P_diag_mc,3));

% Plot mean of estimated states
plot(ax_state(3),t_range,x_est_mean(1,:)','Color',c4,'LineWidth',2,'DisplayName','x1(t) Estimated Mean EnKF');
plot(ax_state(4),t_range,x_est_mean(2,:)','Color',c4,'LineWidth',2,'DisplayName','x2(t) Estimated Mean EnKF');

% Start Monte-Carlo labels
plot(ax_error(1),nan,nan,'Color',c1w,'LineWidth',1,'DisplayName',['x1(t) Error ',mc_str]);
plot(ax_error(2),nan,nan,'Color',c1w,'LineWidth',1,'DisplayName',['x2(t) Error ',mc_str]);
plot(ax_error(3),nan,nan,'Color',c1w,'LineWidth',1,'DisplayName',['Innovation Error ',mc_str]);

% Plot estimation errors for all Monte-Carlo runs
for iter_mc = 1:N_mc
    plot(ax_error(1),t_range,x_err_mc(1,:,iter_mc),'Color',c1w,'LineWidth',1,'HandleVisibility','off');
    plot(ax_error(2),t_range,x_err_mc(2,:,iter_mc),'Color',c1w,'LineWidth',1,'HandleVisibility','off');
    plot(ax_error(3),t_range,innov_mc(:,iter_mc),'Color',c1w,'LineWidth',1,'HandleVisibility','off');
end

% Plot mean of estimation errors
plot(ax_error(1),t_range,x_err_mean(1,:),'Color',c1,'LineWidth',2,'DisplayName','x1(t) Error Mean EnKF');
plot(ax_error(2),t_range,x_err_mean(2,:),'Color',c1,'LineWidth',2,'DisplayName','x2(t) Error Mean EnKF');
plot(ax_error(3),t_range,innov_mean,'Color',c1,'LineWidth',2,'DisplayName','Innovation Mean EnKF');

% Plot theoretical standard deviation
plot(ax_error(1),t_range,P_std(1,:),'--','Color',c0,'LineWidth',2,'DisplayName','sqrt(P1) EnKF');
plot(ax_error(1),t_range,-P_std(1,:),'--','Color',c0,'LineWidth',2,'HandleVisibility','off');
plot(ax_error(2),t_range,P_std(2,:),'--','Color',c0,'LineWidth',2,'DisplayName','sqrt(P2) EnKF');
plot(ax_error(2),t_range,-P_std(2,:),'--','Color',c0,'LineWidth',2,'HandleVisibility','off');

% Plot computed RMS
plot(ax_error(1),t_range,x_est_rms(1,:),'Color',c4,'LineWidth',2,'DisplayName','x1(t) RMS EnKF');
plot(ax_error(1),t_range,-x_est_rms(1,:),'Color',c4,'LineWidth',2,'HandleVisibility','off');
plot(ax_error(2),t_range,x_est_rms(2,:),'Color',c4,'LineWidth',2,'DisplayName','x2(t) RMS EnKF');
plot(ax_error(2),t_range,-x_est_rms(2,:),'Color',c4,'LineWidth',2,'HandleVisibility','off');

% Link plots axes
linkaxes(ax_state,'x');
linkaxes(ax_error,'xy');

%% Compare CDEnKF with CDUKF and CDEFK

% Store EnKF statistics
x_err_mean_EnKF = x_err_mean;
innov_mean_EnKF = innov_mean;
x_est_rms_EnKF  = x_est_rms;
P_std_EnKF      = P_std;

% Run EKF
tic;
for iter_mc = 1:N_mc
    [x_est, P, innov] = extended_kalman_filter(f,F,h,H,P_dot,G,R,Q,u_mc(:,iter_mc),y_mc(:,iter_mc),t_step,x_bar,P_bar);
    x_est_mc(:,:,iter_mc) = x_est;
    innov_mc(:,iter_mc) = innov;
    P_diag_mc(:,:,iter_mc) = [P(1,1,:); P(2,2,:)];
end
EKF_dt = toc;
disp(['EKF time to run: ',num2str(EKF_dt), ' s']);

% Calculate EKF statistics
x_err_mc   = x_est_mc - x_mc;
x_err_mean = mean(x_err_mc,3);
x_est_mean = mean(x_est_mc,3);
innov_mean = mean(innov_mc,2);
x_est_rms  = sqrt(mean(x_err_mc.^2, 3));
P_std = sqrt(mean(P_diag_mc,3));

% Store EKF statistics
x_err_mean_EKF = x_err_mean;
innov_mean_EKF = innov_mean;
x_est_rms_EKF  = x_est_rms;
P_std_EKF      = P_std;

% Run UKF
tic;
for iter_mc = 1:N_mc
    [x_est, P, innov] = unscented_kalman_filter(f,h,G,R,Q,u_mc(:,iter_mc),y_mc(:,iter_mc),t_step,x_bar,P_bar);
    x_est_mc(:,:,iter_mc) = x_est;
    innov_mc(:,iter_mc) = innov;
    P_diag_mc(:,:,iter_mc) = [P(1,1,:); P(2,2,:)];
end
UKF_dt = toc;
disp(['EKF time to run: ',num2str(UKF_dt), ' s']);

% Calculate UKF statistics
x_err_mc   = x_est_mc - x_mc;
x_err_mean = mean(x_err_mc,3);
x_est_mean = mean(x_est_mc,3);
innov_mean = mean(innov_mc,2);
x_est_rms  = sqrt(mean(x_err_mc.^2, 3));
P_std = sqrt(mean(P_diag_mc,3));

% Store UKF statistics
x_err_mean_UKF = x_err_mean;
innov_mean_UKF = innov_mean;
x_est_rms_UKF  = x_est_rms;
P_std_UKF      = P_std;

% Plot differences
plot(ukf_vs_ekf(1),t_range,x_err_mean_EKF(1,:),'Color',c2,'LineWidth',2,'DisplayName','x1(t) Error EKF');
plot(ukf_vs_ekf(1),t_range,x_err_mean_UKF(1,:),'Color',c3,'LineWidth',2,'DisplayName','x1(t) Error UKF');
plot(ukf_vs_ekf(1),t_range,x_err_mean_EnKF(1,:),'Color',c4,'LineWidth',2,'DisplayName','x1(t) Error EnKF');
% plot(ukf_vs_ekf(1),t_range,P_std_EKF(1,:),'--','Color',c2c4,'LineWidth',2,'DisplayName','sqrt(P1) EKF');
% plot(ukf_vs_ekf(1),t_range,P_std_UKF(1,:),'--','Color',c1c3,'LineWidth',2,'DisplayName','sqrt(P1) UKF');
% plot(ukf_vs_ekf(1),t_range,x_est_rms_EKF(1,:),'Color',c2c4,'LineWidth',2,'DisplayName','x1(t) RMS EKF');
% plot(ukf_vs_ekf(1),t_range,x_est_rms_UKF(1,:),'Color',c1c3,'LineWidth',2,'DisplayName','x1(t) RMS UKF');

plot(ukf_vs_ekf(2),t_range,x_err_mean_EKF(2,:),'Color',c2,'LineWidth',2,'DisplayName','x2(t) Error EKF');
plot(ukf_vs_ekf(2),t_range,x_err_mean_UKF(2,:),'Color',c3,'LineWidth',2,'DisplayName','x2(t) Error UKF');
plot(ukf_vs_ekf(2),t_range,x_err_mean_EnKF(2,:),'Color',c4,'LineWidth',2,'DisplayName','x1(t) Error EnKF');
% plot(ukf_vs_ekf(2),t_range,P_std_EKF(2,:),'--','Color',c2c4,'LineWidth',2,'DisplayName','sqrt(P1) EKF');
% plot(ukf_vs_ekf(2),t_range,P_std_UKF(2,:),'--','Color',c1c3,'LineWidth',2,'DisplayName','sqrt(P1) UKF');
% plot(ukf_vs_ekf(2),t_range,x_est_rms_EKF(2,:),'Color',c2c4,'LineWidth',2,'DisplayName','x1(t) RMS EKF');
% plot(ukf_vs_ekf(2),t_range,x_est_rms_UKF(2,:),'Color',c1c3,'LineWidth',2,'DisplayName','x1(t) RMS UKF');

plot(ukf_vs_ekf(3),t_range,innov_mean_EKF,'Color',c2,'LineWidth',2,'DisplayName','Innovation EKF');
plot(ukf_vs_ekf(3),t_range,innov_mean_UKF,'Color',c3,'LineWidth',2,'DisplayName','Innovation UKF');
plot(ukf_vs_ekf(3),t_range,innov_mean_EnKF,'Color',c4,'LineWidth',2,'DisplayName','Innovation EnKF');

%% Continuous-Discrete Ensemble Kalman Filter Definition
function [x_est, P, innov] = ensemble_kalman_filter(f,h,G,R,Q,u,y,dt,x11,P11,N_samples)
    % Matrix sizes
    N_states = length(x11);
    N_steps  = length(y)-1;
    N_mes    = size(R, 1);

    % Initialization (k=1)
    xkk = x11;
    Pkk = P11;
    Xkk = xkk + sqrt(Pkk)*randn(N_states, N_samples);


    % Lists to store estimative and filter covariance
    x_est      = zeros(N_states,N_steps+1);
    P          = zeros(N_states,N_states,N_steps+1);
    innov      = zeros(1,N_steps+1);
    x_est(:,1) = xkk;
    P(:,:,1)   = Pkk;

    for k = 1:N_steps
        % Generate new random samples for forecasting
        W = sqrt(Q)*randn(N_states, N_samples);
        V = sqrt(R)*randn(N_mes, N_samples);

        % Forcast step
        Xk1k  = zeros(N_states,N_samples);
        Yk1k  = zeros(N_mes,N_samples);
        for i = 1:N_samples
            Xk1k(:,i) = rk4_x(f,Xkk(:,i),u(k),G*W(:,i),dt);
            Yk1k(:,i) = h(Xk1k(:,i)) + V(:,i);
        end

        % Prediction
        xk1k   = mean(Xk1k,2);
        dx     = Xk1k - xk1k;
        Pk1k   = (dx*dx')/(N_samples-1); % Unused, kept here for documentation
        yk1k   = mean(Yk1k,2);
        dy     = Yk1k - yk1k;
        Pyk1k  = (dy*dy')/(N_samples-1);
        Pxyk1k = (dx*dy')/(N_samples-1);

        % Update
        Kk1   = Pxyk1k / Pyk1k;
        Xk1k1 = Xk1k + Kk1 * (y(k+1)-Yk1k);
        xk1k1 = mean(Xk1k1,2);
        dx    = Xk1k1 - xk1k1;
        Pk1k1 = (dx*dx')/(N_samples-1);

        % Refresh variables
        Xkk = Xk1k1;

        % Save results
        x_est(:,k+1) = xk1k1;
        P(:,:,k+1)   = Pk1k1;
        innov(k+1)   = y(k+1) - mean(Yk1k);
    end
end

%% Continuous-Discrete Unscented Kalman Filter Definition
function [x_est, P, innov] = unscented_kalman_filter(f,h,G,R,Q,u,y,dt,x11,P11)
    % Initialization (k=1)
    xkk = x11;
    Pkk = P11;

    % Matrix sizes
    N_states = length(x11);
    N_steps  = length(y)-1;
    N_mes    = size(R, 1);
    N_a      = 2*N_states + N_mes;

    % Lists to store estimative and filter covariance
    x_est      = zeros(N_states,N_steps+1);
    P          = zeros(N_states,N_states,N_steps+1);
    innov      = zeros(1,N_steps+1);
    x_est(:,1) = xkk;
    P(:,:,1)   = Pkk;

    % UKF parameters (must be defined assuming the augmented state)
    kappa  = 3 - N_a; % scale parameter
    gamma2 = N_a + kappa;
    gamma  = sqrt(gamma2);
    rho    = [kappa/gamma2, repmat(0.5/gamma2, 1, 2*N_a)];

    for k = 1:N_steps
        % Augmented state
        ax_bar_k = [xkk; zeros(N_states,1); zeros(N_mes,1)];
        aPk      = blkdiag(Pkk,Q,R) + 1e-10*eye(N_a);

        % Generate augmented sigma-points
        aux_sqrt = gamma * chol(aPk, 'lower');
        aX_sigma = [ax_bar_k, ax_bar_k + aux_sqrt, ax_bar_k - aux_sqrt];

        % Get individual sigma-points
        X_sigma = aX_sigma(1            : N_states ,:);
        W_sigma = aX_sigma(N_states+1   : 2*N_states ,:);
        V_sigma = aX_sigma(2*N_states+1 : end ,:);

        % Propagate the sigma points through integration
        X_sigma_k1k = zeros(N_states, 2*N_a+1);
        Y_sigma_k1k = zeros(N_mes, 2*N_a+1);
        for i = 1:2*N_a+1
            X_sigma_k1k(:,i) = rk4_x(f,X_sigma(:,i),u(k),G*W_sigma(:,i),dt);
            Y_sigma_k1k(:,i) = h(X_sigma_k1k(:,i)) + V_sigma(:,i);
        end

        % Prediction
        xk1k = 0;
        yk1k = 0;
        for i = 1:2*N_a+1
            xk1k = xk1k + rho(i)*X_sigma_k1k(:,i);
            yk1k = yk1k + rho(i)*Y_sigma_k1k(:,i);
        end

        innovk = y(k+1) - yk1k;

        Pk1k   = zeros(N_states);
        Pyk1k  = zeros(N_mes);
        Pxyk1k = zeros(N_states,N_mes);
        for i = 1:2*N_a+1
            dx = X_sigma_k1k(:,i)-xk1k;
            dy = Y_sigma_k1k(:,i)-yk1k;
            Pk1k   = Pk1k   + rho(i)*(dx*dx');
            Pyk1k  = Pyk1k  + rho(i)*(dy*dy');
            Pxyk1k = Pxyk1k + rho(i)*(dx*dy');
        end
        Pyk1k = Pyk1k + 1e-10*eye(N_mes);

         % Update
        Kk1   = Pxyk1k / Pyk1k;
        xk1k1 = xk1k + Kk1 * innovk;
        Pk1k1 = Pk1k - Pxyk1k/Pyk1k * Pxyk1k';
        Pk1k1 = 0.5 * (Pk1k1 + Pk1k1'); % make it symetric

        % Refresh variables
        xkk = xk1k1;
        Pkk = Pk1k1;

        % Save results
        x_est(:,k+1) = xk1k1;
        P(:,:,k+1)   = Pk1k1;
        innov(k+1)   = innovk;
    end
end

%% Continuous-Discrete Extended Kalman Filter Definition
function [x_est, P, innov] = extended_kalman_filter(f,F,h,H,P_dot,G,R,Q,u,y,dt,x11,P11)
    % Initialization (k=1)
    xkk = x11;
    Pkk = P11;

    % Matrix sizes
    N_states = length(x11);
    N_steps  = length(y)-1;

    % Lists to store estimative and filter covariance
    x_est      = zeros(N_states,N_steps+1);
    P          = zeros(N_states,N_states,N_steps+1);
    innov      = zeros(1,N_steps+1);
    x_est(:,1) = xkk;
    P(:,:,1)   = Pkk;

    for k = 1:N_steps
        % Integration and Prediction
        xk1k   = rk4_x(f,xkk,u(k),0,dt);
        Pk1k   = rk4_P(P_dot,xk1k,F,Pkk,G,Q,dt);
        yk1k   = h(xk1k);
        innovk = y(k+1) - yk1k;
        Pyk1k  = H * Pk1k * H' + R;
        Pxyk1k = Pk1k * H';

        % Update
        Kk1   = Pxyk1k / Pyk1k;
        xk1k1 = xk1k + Kk1 * innovk;
        Pk1k1 = Pk1k - Pxyk1k/Pyk1k * Pxyk1k';

        % Refresh variables
        xkk = xk1k1;
        Pkk = Pk1k1;

        % Save results
        x_est(:,k+1) = xk1k1;
        P(:,:,k+1)   = Pk1k1;
        innov(k+1)   = innovk;
    end
end

%% Runge-Kutta 4 for the state x
function xk1 = rk4_x(f,xk,uk,wk,dt)
    x1 = xk;
    k1 = f(x1, uk) + wk;
    x2 = x1 + 0.5 * dt * k1;
    k2 = f(x2, uk) + wk;
    x3 = x1 + 0.5 * dt * k2;
    k3 = f(x3, uk) + wk;
    x4 = x1 + dt * k3;
    k4 = f(x4, uk) + wk;

    xk1 = x1 + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
end

%% Runge-Kutta 4 for the covariance P
function Pk1 = rk4_P(f,xk,F,P,G,Q,dt)
    P1 = P;
    k1 = f(xk,F,P1,G,Q);
    P2 = P1 + 0.5 * dt * k1;
    k2 = f(xk,F,P2,G,Q);
    P3 = P1 + 0.5 * dt * k2;
    k3 = f(xk,F,P3,G,Q);
    P4 = P1 + dt * k3;
    k4 = f(xk,F,P4,G,Q);

    Pk1 = P1 + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
end

%% Plot Setup just to improve readability
function ax = plot_setup(title_str,x_str,y_str)
    figure;

    n_plots = length(title_str);
    ax = zeros(n_plots,1);
    for i = 1:n_plots
        ax(i) = subplot(n_plots,1,i);
        hold on;
        grid on;
        legend;
        title(title_str{i});
        xlabel(x_str{i});
        ylabel(y_str{i});
    end
end
