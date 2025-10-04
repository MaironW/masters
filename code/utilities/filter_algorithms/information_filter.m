% Simple Discrete Information Filter implementation
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