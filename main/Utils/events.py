# Change parameter values over time

import numpy as np
from copy import deepcopy

# Default eeents definition (time, value)
# The simulation will hold the lastest value based on the current time step
events_sequence = {
    "DYN.DYN_ATT.BOFq_SSB" : [
        (0, np.array([1,0,0,0])),
    ],
    "DYN.DYN_ATT.SCw_BOF" : [
        (0, np.deg2rad([0,0,0])),
    ],
    "SEN.SEN_STR.STRenableflg" : [
        (100,   1),
        (6000,  1),
        (36000, 1),
    ],
    "SEN.SEN_PSR.PSRenableflg" : [
        (0,     1),
        (3000,  1),
        (36000, 1),
    ],
    "SEN.SEN_CMB.CMBenableflg" : [
        (0,     1),
        (3000,  1),
        (36000, 1),
    ],
    "NAV.NAV_CEL.NAV_CELenableflg" : [
        (0, 1),
    ],
    "NAV.NAV_PSR.NAV_PSRenableflg" : [
        (0, 1),
    ],
    "NAV.NAV_CMB.NAV_CMBenableflg" : [
        (0, 1),
    ]
}

# Build a full per-time-step dictionary of parameter values
def build_events_table(t_start, t_end, dt, events=events_sequence):
    # Intialize the value state for all parameters
    current_values = {}
    for var, changes in events.items():
        # Sort each variable's events by time
        sorted_changes = sorted(changes, key=lambda x: x[0])
        current_values[var] = sorted_changes[0][1] # initial value
    # Precompute event times for faster lookup
    event_times = {
        var: [t for t, _ in sorted(changes, key=lambda x: x[0])]
        for var, changes in events.items()
    }
    event_values = {
        var: [v for _, v in sorted(changes, key=lambda x: x[0])]
        for var, changes in events.items()
    }
    # Create events_table
    events_table = {}
    # Iterate over simulation time
    for t in range(t_start, t_end + 1, dt):
        # Update current values if any event occurs
        for var, times in event_times.items():
            for i, event_t in enumerate(times):
                if event_t <= t:
                    current_values[var] = event_values[var][i]
        # Group values by module
        grouped = {}
        for full_key, value in current_values.items():
            parts = full_key.split(".")
            if len(parts) < 2:
                continue
            level1 = parts[0] # DYN, SEN, NAV
            level2 = parts[1] # DYN.###
            level3 = parts[2] # DYN.###.param
            if level1 not in grouped:
                grouped[level1] = {}
            if level2 not in grouped[level1]:
                grouped[level1][level2] = {}
            grouped[level1][level2][level3] = value
        events_table[t] = deepcopy(grouped)
    return events_table
