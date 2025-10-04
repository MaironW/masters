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