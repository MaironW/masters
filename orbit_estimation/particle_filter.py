from numba import njit

import numpy as np
from scipy.stats import multivariate_normal

# Adapted from aleskandarhabe.com
@njit
def systematic_resampling(weight_array):
    # N is the total number of samples
    N = len(weight_array)
    divide_by_N = 1/N # to make the code a bit faster
    # cummulative sum of weights
    c_values = np.cumsum(weight_array)
    # starting random point
    starting_point = np.random.uniform(low=0, high=divide_by_N)

    # this list stores indices of resampled states
    resampled_index=[0]*N
    for j in range(N):
        current_point = starting_point+j*divide_by_N
        s=0
        # get the index of the first point in the cumulative sum lower than the current point
        # somehow the method below is faster than numpy and bisect
        while (current_point > c_values[s]):
            s=s+1
        resampled_index[j] = s
    return resampled_index

def particle_filter(F, B, Z, u, R, Q, dt, n_particles):

    # Matrix dimensions
    n_states = Q.shape[0] # state dimension
    n_iters  = Z.shape[0] # number of steps

    # Create distribuition for the process
    # allow_singular=True must be used to allow ill_conditioned matrix
    # TODO: learn how to make the matrix Q not ill conditioned
    process_distribuition = multivariate_normal(mean=[0]*n_states, cov=Q, allow_singular=True)

    # Lists to store the states and weights
    state_list     = [0]*n_iters
    state_list[0]  = Z[0] # first state equal to the initial measurement
    weight_list    = [0]*n_iters
    weight_list[0] = (1/n_particles)*np.ones((1, n_particles))

    # Allocate result matrices
    state_list      = np.zeros((n_iters, n_states, n_particles))
    mean_state_list = np.zeros((n_iters, n_states))
    weight_list     = np.zeros((n_iters, n_particles))
    new_weights     = np.zeros((1, n_particles))

    # Set initial state and weight values for each particle
    state_init  = np.zeros((n_states,n_particles))
    for i in range(n_states):
        state_range = [Z[0, i]-R[i,i], Z[0, i]+R[i,i]]
        state_init[i,:] = np.random.uniform(state_range[0], state_range[1], size=n_particles)

    states  = state_init
    weights = (1/n_particles)*np.ones((1, n_particles))

    # Save the initial values
    state_list[0]  = states
    weight_list[0] = weights

    # Run filter
    for i in range(n_iters):
        print(i)

        # STEP 1: Predict
        # Apply the state transition probability to the particles
        predicted_states = F(states, dt)
        process_random_variables = process_distribuition.rvs(size=n_particles).transpose()
        new_states = predicted_states + np.matmul(B,u[i]).reshape(n_states,1) + process_random_variables

        # STEP 2: Update the weights
        # [MODIFICATION] only run when a new measurement arrives
        # If no measurement is available, do not update the weights
        if not np.array_equal(Z[i],Z[i-1]):
            for j in range(n_particles):
                # Create a normal distribuition on the basis of the current particle state
                # Measurement distribuition
                measurement_distribuition = multivariate_normal(mean=new_states[:, j], cov=R)
                new_weights[:, j] = measurement_distribuition.pdf(Z[i])*weights[:, j]
            new_weights_std = new_weights/(new_weights.sum())
            # STEP 3: Resample
            tmp = [val**2 for val in new_weights_std]
            Neff = 1/np.array(tmp).sum()

            if Neff < (n_particles//3):
                # resampled_state_index = np.random.choice(np.arange(n_particles), n_particles, p=new_weights_std[0])
                resampled_state_index = systematic_resampling(new_weights_std[0,:])
                new_states      = new_states[:, resampled_state_index]
                new_weights_std = (1/n_particles)*np.ones((1, n_particles))

        # compute the mean state
        # for the next iteration
        states  = new_states
        weights = new_weights_std

        # save the states and weights
        mean_state_list[i] = (weights*states).sum(axis=1)
        state_list[i]      = states
        weight_list[i]     = weights

    return mean_state_list, state_list, weight_list