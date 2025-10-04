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