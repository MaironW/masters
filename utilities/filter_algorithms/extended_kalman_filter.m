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