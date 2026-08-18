# Sample Monte Carlo script
# Run from repository root with:
#   python3 -m Tests.test_mc

import os
import json
import time
import numpy as np
import multiprocessing as mp

from PPC   import PPC
from Utils import simulation

# Monte Carlo function for multiprocessing
def run_mc_case(args):
    i, output_dir = args
    print(f"[MC] Run {i+1}/{mc_steps}")
    run_id  = f"{i:04d}"
    run_dir = os.path.join(output_dir, f"run_{run_id}")
    os.makedirs(run_dir, exist_ok=True)

    # Randomization
    seed = np.random.randint(0, 1e9)
    np.random.seed(seed)

    # Define parameters distribuition
    P_ini       = np.diag(np.random.uniform(0.5, 2.0, size=6))
    pos_ini_err = np.random.normal(0, 10,   size=3)
    vel_ini_err = np.random.normal(0, 0.01, size=3)
    x_ini_err   = np.hstack([pos_ini_err, vel_ini_err])
    x_ini_true  = np.array([1.62228612e+08, 2.42235379e+07, 1.78903662e+07, -0.32883677, 27.86282899, 13.43536119])
    x_est_ini   = x_ini_true + x_ini_err

    # Parameter override
    par_override = {
        "NAV": {
            "NAV_EKF": {
                "P_ini"     : P_ini,
                "x_est_ini" : x_est_ini,
            }
        }
    }

    # Simulation override
    SIM_par_override = {
        "log_save_path": os.path.join(run_dir, "log"),
        "log_load_path": os.path.join(output_dir, "baseline"),
        "log_save": True,
        "log_load": ["DYN"]
    }

    # Save parameters
    metadata = {
        "seed"         : seed,
        "par_override" : par_override,
    }

    with open(os.path.join(run_dir, "params.json"), "w") as f:
        json.dump(PPC.serialize_json(metadata), f, indent=4)

    # Run simulation
    simulation.run(
        SIM_par_override = SIM_par_override,
        par_override     = par_override,
        run_id           = run_id,
        log_time         = False,
    )

# Main
if __name__ == "__main__":

    script_start_time = time.perf_counter()

    mc_steps   = 2
    output_dir = "Logs/test_mc"
    os.makedirs(output_dir, exist_ok=True)

    # Run baseline simulation (To avoid running DYN every execution)
    print(f"[MC] Baseline Run")
    SIM_par_override = {
        "log_save_path": output_dir+"/baseline",
        "log_save": True,
        "log_load": []
    }
    simulation.run(SIM_par_override=SIM_par_override, run_id="")

    # Parallel execution
    n_proc = mp.cpu_count()

    with mp.Pool(processes=n_proc) as pool:
        pool.map(run_mc_case, [(i, output_dir) for i in range(mc_steps)])

    # Stop simulation
    simulation.stop()

    # Log time
    script_end_time = time.perf_counter()
    runtime = script_end_time - script_start_time
    print(f"[MC] Total Run Time: {runtime:.3f} s")
