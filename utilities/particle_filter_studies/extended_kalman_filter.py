import numpy as np
from numpy.linalg import inv

def extended_kalman_filter(JF,B,JH,Z,u,R,Q):

    # matrix dimensions
    nx = Q.shape[0]
    ny = R.shape[0]
    nt = Z.shape[0]

    # allocate identity matrix for re-use
    I = np.eye(nx)

    # initial estimative
    x_init = Z[0]
    P_init = 0.1 * np.eye(len(x_init))  # small initial prediction error

    # allocate result matrices
    x_pred = np.zeros((nt, nx))      # prediction of state vector
    P_pred = np.zeros((nt, nx, nx))  # prediction error covariance matrix
    x_est = np.zeros((nt, nx))       # estimation of state vector
    P_est = np.zeros((nt, nx, nx))   # estimation error covariance matrix
    K = np.zeros((nt, nx, ny))       # Kalman Gain

    # set initial prediction
    x_pred[0] = x_init
    P_pred[0] = P_init
    x_est[0] = x_init
    P_est[0] = P_init

    # for each time-step...
    for i in range(1,nt):

        # prediction stage
        x_pred[i] = np.matmul(JF,x_est[i-1]) + np.matmul(B,u[i])
        P_pred[i] = np.matmul(np.matmul(JF,P_est[i-1]),JF.T) + Q

        # estimation stage
        # [MODIFICATION] only run when a new measurement arrives
        if not np.array_equal(Z[i],Z[i-1]):
            y = Z[i] - np.matmul(JH,x_pred[i])

            # Kalman gain
            S = np.matmul(np.matmul(JH,P_pred[i]),JH.T) + R
            K[i] = np.matmul(np.matmul(P_pred[i],JH.T), inv(S))

            # prediction
            x_est[i] = x_pred[i] + np.matmul(K[i],y)
            P_est[i] = np.matmul((I - np.matmul(K[i],JH)),P_pred[i])
            # P_est[i] = np.matmul(np.matmul(I-K[i],P_pred[i]),(I-K[i]).T) + np.matmul(np.matmul(K[i],R),K[i].T)
        else:
            x_est[i] = x_pred[i]
            P_est[i] = P_pred[i]

    return x_est