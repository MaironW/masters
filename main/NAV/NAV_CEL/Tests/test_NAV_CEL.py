# NAV_CEL Module Test
# Generates the sensor field of view when the spacecraft is aways pointing to Mars
# Compares the measurements with the simulated stars from NAV_CEL
# Run from repository root with:
#   python3 -m NAV.NAV_CEL.Tests.test_NAV_CEL

import numpy as np
from DYN.DYN         import DYN
from SEN.SEN         import SEN
from NAV.NAV         import NAV
from PPC             import PPC
from PPC.PPC_plots   import PPC_plots
from Utils           import spice
from Utils           import events
from Utils           import integrator
from Utils           import quaternions
from Utils.constants import CONSTANTS_par

#########
# SETUP #
#########

SIM_par = {
    # Time parameters
    "dt"           : 600,       # [s] 10 min
    "time_start"   : 0,         # [s]
    "time_end"     : 3600*24*5, # [s] 5 days

    # Log parameters
    "DYN_log_save" : False,
    "DYN_log_load" : False,
    "DYN_log_path" : "Logs/DYN",

    # List of default plots
    "PPC_plot_list" : [
        "SEN_STR",
        "NAV_STR",
        "NAV_CEL",
    ],
}

# Load data or run new simulation
DYN_log_save = SIM_par["DYN_log_save"]
DYN_log_load = SIM_par["DYN_log_load"]
DYN_log_path = SIM_par["DYN_log_path"]

# Simulation parameters
sim_dt         = SIM_par["dt"]
sim_time_start = SIM_par["time_start"]
sim_time_end   = SIM_par["time_end"]
sim_time       = sim_time_start
step           = 0
n_steps        = int((sim_time_end - sim_time_start)/sim_dt) + 1

# Set events
events_sequence = events.events_sequence
events_sequence["SEN.SEN_STR.STRenableflg"] = [
        (100,   0),
        (6000,  1),
        (12000, 0),
        (13000, 1),
]

# Initialize inputs table
events_table = events.build_events_table(sim_time_start, sim_time_end, sim_dt, events=events_sequence)

##################
# MRO PARAMETERS #
##################

# Load NASA's MRO data with SPICE
# https://naif.jpl.nasa.gov/pub/naif/pds/data/mro-m-spice-6-v1.0/mrosp_1000/data
kernel_dir = "Utils/kernels/"
spice.load_kernel(kernel_dir + "mro_cruise.bsp")
spice.load_kernel(kernel_dir + "mro_sclkscet_00021_65536.tsc")

# Get Spacecraft initial condition
# Set time so that spacecraft starts closer to Mars
time_TDB_ini = spice.get_time("2006-03-05 T00:00:00")

MROpos_SSB_ini, MROvel_SSB_ini = spice.get_state("MRO", time_TDB_ini)

######################
# INITIALIZE MODULES #
######################

DYN_TIME_par = {
    "time_SIM_ini" : 0,
    "time_TDB_ini" : time_TDB_ini,
}

DYN_TRA_par = {
    "ref_elements"  : "rvi",
    "BODY_ini"      : "SUN",
    "sma_ini"       : None,
    "ecc_ini"       : None,
    "incl_ini"      : None,
    "raan_ini"      : None,
    "argp_ini"      : None,
    "tano_ini"      : None,
    "SCpos_SSB_ini" : MROpos_SSB_ini, # [km]   Initial position relative to SSB frame
    "SCvel_SSB_ini" : MROvel_SSB_ini, # [km/s] Initial velocity relative to SSB frame
}

field_of_view     = 30*CONSTANTS_par["deg2rad_cst"]
cos_field_of_view = np.cos(field_of_view)

SEN_STR_par = {
    "field_of_view"     : field_of_view,
    "cos_field_of_view" : cos_field_of_view,
}

par_override = {
    "DYN_TIME" : DYN_TIME_par,
    "DYN_TRA"  : DYN_TRA_par,
    "SEN_STR"  : SEN_STR_par,
}

# Initialize DYN modules with overrided time
DYN_obj = DYN(par_override)

# Set only modules needed for the test
DYN_obj.modules = [
    DYN_obj.DYN_TIME,
    DYN_obj.DYN_EPH,
    DYN_obj.DYN_TRA,
    DYN_obj.DYN_GRV,
    DYN_obj.DYN_ATT,
    DYN_obj.DYN_STR,
    DYN_obj.DYN_PSR,
]

# Initialize SEN modules with parameter override
SEN_obj = SEN(DYN_obj, par_override)

# Normally initialize NAV modules
NAV_obj = NAV(SEN_obj)

# Set only modules needed for the test
NAV_obj.modules = [
    NAV_obj.NAV_EPH,
    NAV_obj.NAV_STR,
    NAV_obj.NAV_CEL,
]

# Initialize timeline
timeline = {
    "DYN" : PPC.init_timeline(DYN_obj.snapshot(), n_steps),
    "SEN" : PPC.init_timeline(SEN_obj.snapshot(), n_steps),
    "NAV" : PPC.init_timeline(NAV_obj.snapshot(), n_steps)
}

########
# LOOP #
########

# Main loop
for step in range(1, n_steps):

    # Load external inputs
    inputs = events_table[sim_time]

    ###################
    # UPDATE ATTITUDE #
    ###################

    # Spacecraft is aways pointing +STRz to Mars
    DYN_last_state = DYN_obj.snapshot()
    MARSpos_SSB = DYN_last_state["DYN_EPH"]["MARSpos_SSB"]
    SCpos_SSB   = DYN_last_state["DYN_TRA"]["SCpos_SSB"]
    STRy_STR    = np.array([0,1,0])
    STRz_STR    = np.array([0,0,1])
    STRq_BOF    = SEN_obj.SEN_STR.par["STRq_BOF"]

    BOFq_STR    = quaternions.qtrans(STRq_BOF)

    # Line of sight from SC to Mars
    MARSpos_SC = MARSpos_SSB - SCpos_SSB
    MARSdir_SC = SEN_obj.SEN_STR.dir_from_pos(MARSpos_SC)

    # STR boresight in the BOF frame
    STRz_BOF = quaternions.qvecprod(BOFq_STR, STRz_STR)

    # Get the quaternion which points STRz_BOF to MARSu_SSB
    BOFq_SSB_tgt = quaternions.vecqvec(STRz_BOF, MARSdir_SC)

    inputs["DYN"]["DYN_ATT"]["BOFq_SSB"] = BOFq_SSB_tgt

    ###################

    # Update algebraic modules first
    DYN_obj.update_algebraic(sim_time, None, inputs)
    # Integrate all dynamic states together
    integrator.rk4_step(sim_time, sim_dt, DYN_obj, inputs)
    # Update SEN
    SEN_obj.update_algebraic(sim_time, DYN_obj, inputs)
    # Update NAV
    NAV_obj.update_algebraic(sim_time, SEN_obj, inputs)

    # Save results into the timeline
    output = {
        "DYN" : DYN_obj.snapshot(),
        "SEN" : SEN_obj.snapshot(),
        "NAV" : NAV_obj.snapshot(),
    }
    PPC.update_timeline(timeline, output, step)

    # Update time
    sim_time += sim_dt

############
# TEARDOWN #
############

# Clear Kernels from memory
spice.clear_kernels()

###########
# RESULTS #
###########

for key in SIM_par["PPC_plot_list"]:
    if key in PPC_plots:
        PPC_plots[key](timeline, DYN_obj, SEN_obj, NAV_obj)
PPC.show_plot()
