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