% MP208
% Computational Exercise 2
% Mairon de Souza Wolniewicz
% Written in MATLAB R2018b

%% Setup
clc; clear; close all;

% Simulation parameters
N_mc    = 100;   % Number of realizations
t_start =  0.0; % [s] Start time
t_end   = 20.0; % [s] End time
t_step  =  0.1; % [s] Sample time (T from the problem statement)
t_range = t_start:t_step:t_end;
N_steps = length(t_range); % Number of steps

% Model parameters
Q = diag([1e-2, 4e-2]);     % Process covariance matrix
R = 1e-2;                   % Measurement variance

x_bar = [1; 0];             % Initial state mean
P_bar = diag([1e-4, 1e-8]); % Initial state covariance

% u parameters
y_bar = 5;
e1 = [1; 0];
e2 = [0; 1];

% State, Control and Measurement Matrices
A = [1 0.1; 0 1];
B = [0.005; 0.1];
C = [1 0];
G = eye(2);

% Preset plots so lines are generated in computation order
ax_y       = plot_setup('System Output','k','y');
ax_x1      = plot_setup('State x1','k','x1');
ax_x2      = plot_setup('State x2','k','x2');
ax_kf_err1 = plot_setup('Kalman Filter x1 Error','k','x1\_err');
ax_kf_est1 = plot_setup('Kalman Filter x1 Estimation','k','x1\_est');
ax_kf_err2 = plot_setup('Kalman Filter x2 Error','k','x2\_err');
ax_kf_est2 = plot_setup('Kalman Filter x2 Estimation','k','x2\_est');
ax_kf_inno = plot_setup('Kalman Filter Innovation','k','y - C*x\_est');
ax_if_err1 = plot_setup('Information Filter x1 Error','k','x1\_err');
ax_if_est1 = plot_setup('Information Filter x1 Estimation','k','x1\_est');
ax_if_err2 = plot_setup('Information Filter x2 Error','k','x2\_err');
ax_if_est2 = plot_setup('Information Filter x2 Estimation','k','x2\_est');
ax_time    = plot_setup('Filter Execution Time','Monte Carlo runs','time (s)');

%% Run Monter-Carlo realizations
y_mc      = zeros(N_steps+1, N_mc);
u_mc      = zeros(N_steps, N_mc);
x_mc      = zeros(2, N_steps+1, N_mc);
x_est_mc  = zeros(2, N_steps+1, N_mc);
P_std     = zeros(2, N_steps+1); % Although stored multiple times, P is not dependent on the measurements

%% a) Simulate the system in the period from 0 to 20 s
% Obtain 10 realizations of {yk} and plot together with y_bar
for iter_mc = 1:N_mc
    x = zeros(2,N_steps+1);
    y = zeros(1,N_steps+1);
    y(1) = 1;
    u = zeros(1,N_steps);

    % Realizations
    w      = sqrt(Q) * randn(2,N_steps);
    v      = sqrt(R) * randn(1,N_steps);
    x(:,1) = x_bar + sqrt(P_bar) * randn(2,1);

    for k = 1:N_steps
        u(k)     = 10*(y_bar - e1'*x(:,k)) - 2*e2'*x(:,k);
        x(:,k+1) = A*x(:,k) + B*u(k) + w(:,k);
        y(k+1)   = C*x(:,k+1) + v(k);
    end
    y_mc(:,iter_mc)   = y;
    u_mc(:,iter_mc)   = u;
    x_mc(:,:,iter_mc) = x;

    plot(ax_y,y,'DisplayName',['y (',num2str(iter_mc),')']);
    plot(ax_x1,x(1,:)','DisplayName',['real x1 (',num2str(iter_mc),')']);
    plot(ax_x2,x(2,:)','DisplayName',['real x2 (',num2str(iter_mc),')']);
end
plot(ax_y,y_bar*ones(1,N_steps+1)','--','LineWidth',2,'DisplayName','y_{bar}');
x_mean = mean(x_mc,3);
plot(ax_x1,x_mean(1,:)','k','LineWidth',2,'DisplayName','x1\_real\_mean');
plot(ax_x2,x_mean(2,:)','k','LineWidth',2,'DisplayName','x2\_real\_mean');

%% b) Simulate a Kalman Filter for all 10 realizations of {yk}
disp('Running Kalman filter')
kf_time = zeros(1,N_mc);
for iter_mc = 1:N_mc
    tic;
    [x_est, P, innov] = kalman_filter(A,B,C,G,R,Q,u_mc(:,iter_mc),y_mc(:,iter_mc),x_bar,P_bar);
    kf_time(iter_mc) = toc;
    x_est_mc(:,:,iter_mc) = x_est;
end

% Get statistics
for k=1:N_steps
    P_std(:,k) = sqrt(diag(P(:,:,k)));
end
x_err_mc   = x_est_mc - x_mc;
x_est_mean = mean(x_est_mc,3);
x_est_rms  = sqrt(mean(x_err_mc.^2, 3));

for iter_mc = 1:N_mc
    plot(ax_kf_err1,x_err_mc(1,:,iter_mc)','DisplayName',['x1\_err (',num2str(iter_mc),')']);
    plot(ax_kf_err2,x_err_mc(2,:,iter_mc)','DisplayName',['x2\_err (',num2str(iter_mc),')']);

    plot(ax_kf_est1,x_est_mc(1,:,iter_mc)','DisplayName',['x1\_est (',num2str(iter_mc),')']);
    plot(ax_kf_est2,x_est_mc(2,:,iter_mc)','DisplayName',['x2\_est (',num2str(iter_mc),')']);
end

plot(ax_kf_est1,x_mean(1,:),'k','LineWidth',2,'DisplayName','x1\_real\_mean');
plot(ax_kf_est2,x_mean(2,:),'k','LineWidth',2,'DisplayName','x2\_real\_mean');
plot(ax_kf_est1,x_est_mean(1,:),'b','LineWidth',2,'DisplayName','x1\_est\_mean');
plot(ax_kf_est2,x_est_mean(2,:),'b','LineWidth',2,'DisplayName','x2\_est\_mean');

plot(ax_kf_err1,x_est_rms(1,:),'r','LineWidth',2,'DisplayName','x1\_err\_RMS');
plot(ax_kf_err2,x_est_rms(2,:),'r','LineWidth',2,'DisplayName','x2\_err\_RMS');
plot(ax_kf_err1,P_std(1,:),'m','LineWidth',2,'DisplayName','P1 std');
plot(ax_kf_err2,P_std(2,:),'m','LineWidth',2,'DisplayName','P2 std');

plot(ax_kf_inno,innov,'DisplayName','Innovation');

%% c) Repeat the item b) using the information filter
disp('Running Information filter')
if_time = zeros(1,N_mc);
for iter_mc = 1:N_mc
    tic;
    [x_est, P] = information_filter(A,B,C,G,R,Q,u_mc(:,iter_mc),y_mc(:,iter_mc),x_bar,P_bar);
    if_time(iter_mc) = toc;
    x_est_mc(:,:,iter_mc)  = x_est;
end

% Get statistics
for k=1:N_steps
    P_std(:,k) = sqrt(diag(P(:,:,k)));
end
x_err_mc   = x_est_mc - x_mc;
x_est_mean = mean(x_est_mc,3);
x_est_rms  = sqrt(mean(x_err_mc.^2, 3));

for iter_mc = 1:N_mc
    plot(ax_if_err1,x_err_mc(1,:,iter_mc)','DisplayName',['x1\_err (',num2str(iter_mc),')']);
    plot(ax_if_err2,x_err_mc(2,:,iter_mc)','DisplayName',['x2\_err (',num2str(iter_mc),')']);

    plot(ax_if_est1,x_est_mc(1,:,iter_mc)','DisplayName',['x1\_est (',num2str(iter_mc),')']);
    plot(ax_if_est2,x_est_mc(2,:,iter_mc)','DisplayName',['x2\_est (',num2str(iter_mc),')']);
end

plot(ax_if_est1,x_mean(1,:)','k','LineWidth',2,'DisplayName','x1\_real\_mean');
plot(ax_if_est2,x_mean(2,:)','k','LineWidth',2,'DisplayName','x2\_real\_mean');
plot(ax_if_est1,x_est_mean(1,:)','b','LineWidth',2,'DisplayName','x1\_est\_mean');
plot(ax_if_est2,x_est_mean(2,:)','b','LineWidth',2,'DisplayName','x2\_est\_mean');

plot(ax_if_err1,x_est_rms(1,:)','r','LineWidth',2,'DisplayName','x1\_err\_RMS');
plot(ax_if_err2,x_est_rms(2,:)','r','LineWidth',2,'DisplayName','x2\_err\_RMS');
plot(ax_if_err1,P_std(1,:),'m','LineWidth',2,'DisplayName','P1 std');
plot(ax_if_err2,P_std(2,:),'m','LineWidth',2,'DisplayName','P2 std');

% Compare Kalman and Information filters execution time
plot(ax_time,kf_time,'DisplayName','Kalman Filter');
plot(ax_time,if_time,'DisplayName','Information Filter');

%% Kalman Filter Definition

% Simple Discrete Kalman Filter implementation
% assuming one estimative per measurement
% assuming all measurements were already collected
function [x_est, P, innov] = kalman_filter(A,B,C,G,R,Q,u,y,x11,P11)
    % Initialization (k=1)
    xkk = x11;
    Pkk = P11;

    % Matrix sizes
    N_states = length(x11);
    N_steps  = length(y)-1;

    % Lists to store estimative and filter covariance
    x_est      = zeros(N_states,N_steps+1);
    P          = zeros(N_states,N_states,N_steps+1);
    innov      = zeros(1,N_steps);
    x_est(:,1) = xkk;
    P(:,:,1)   = Pkk;

    for k = 1:N_steps
        % Prediction
        xk1k   = A * xkk + B*u(k);
        yk1k   = C * xk1k;
        innov(k) = y(k+1) - yk1k;
        Pk1k   = A * Pkk * A' + G*Q*G';
        Pyk1k  = C * Pk1k * C' + R;
        Pxyk1k = Pk1k * C';

        % Update
        Kk1 = Pxyk1k / Pyk1k;
        xk1k1 = xk1k + Kk1 * innov(k);
        I = eye(size(Pk1k));
        % Use Joseph formula to avoid errors in the covariance matrix
%         Pk1k1 = Pk1k - Pxyk1k / Pyk1k * Pxyk1k';
        Pk1k1 = (I - Kk1 * C) * Pk1k * (I - Kk1 * C)' + Kk1 * R * Kk1';

        % Refresh variables
        xkk = xk1k1;
        Pkk = Pk1k1;

        % Save results
        x_est(:,k+1) = xk1k1;
        P(:,:,k+1)   = Pk1k1;
    end
end

%% Information Filter Definition

function [x_est, P] = information_filter(A,B,C,G,R,Q,u,y,x11,P11)
    % Matrix sizes
    N_states = length(x11);
    N_steps  = length(y)-1;

    % Matrix initialization
    A_inv = inv(A);
    Q_inv = inv(Q);
    R_inv = inv(R);
    I = eye(N_states);

    % Initialization (k=1)
    xkk = x11;
    Pkk = P11;

    % Get information
    Lkk = inv(Pkk);
    zkk = Lkk*xkk;

    % Lists to store estimative and filter covariance
    x_est      = zeros(N_states,N_steps);
    x_est(:,1) = xkk;
    P          = zeros(N_states,N_states,N_steps);
    P(:,:,1)   = Pkk;

    for k = 1:N_steps
        % Predicton
        PIk = A_inv' * Lkk * A_inv;
        Kk  = PIk * G / (G' * PIk * G + Q_inv);
        aux = (I - Kk * G');
        zk1k = aux * (A_inv' * zkk + PIk * B * u(k));
        Lk1k = aux * PIk;

        % Update
        zk1k1 = zk1k + C' * R_inv * y(k+1);
        Lk1k1 = Lk1k + C' * R_inv * C;

        % Refresh variables
        zkk = zk1k1;
        Lkk = Lk1k1;

        % Save results
        Pk1k1 = inv(Lk1k1);
        xk1k1 = Pk1k1 * zk1k1;

        x_est(:,k+1) = xk1k1;
        P(:,:,k+1)   = Pk1k1;
    end
end

%% Plot Setup just to improve readability
function ax = plot_setup(title_str,x_str,y_str)
    figure;
    ax = axes;
    hold on;
    grid on;
    legend;
    title(title_str);
    xlabel(x_str);
    ylabel(y_str);
end
