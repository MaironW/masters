# DYN_TRA Module Test
# Use the main code solver to propagate spacecraft trajectory
# Compares the result using data from Spice

from ...DYN_TIME  import DYN_TIME
from ...DYN_SUN   import DYN_SUN
from ...DYN_EARTH import DYN_EARTH
from ...DYN_MARS  import DYN_MARS
from ...DYN_GRV   import DYN_GRV
from ...DYN_TRA   import DYN_TRA
from PPC   import PPC
from Utils import spice
from Utils import integrator
from Utils.constants import CONSTANTS_par

# Simulation parameters
sim_dt         = 60*10   # [s] 10 min
sim_time_start = 0       # [s]
sim_time_end   = 3600*24*30 # [s] 30 days
sim_time       = sim_time_start # [s]
step           = 0
n_steps        = int((sim_time_end - sim_time_start)/sim_dt) + 1

# Initialize outputs with correct date/time
DYN_TIME_out  = DYN_TIME.initialize()

# Load NASA's MRO data with SPICE
# https://naif.jpl.nasa.gov/pub/naif/pds/data/mro-m-spice-6-v1.0/mrosp_1000/data
kernel_dir = "DYN/DYN_TRA/Tests/kernels/"
spice.load_kernel(kernel_dir + "mro_cruise.bsp")
spice.load_kernel(kernel_dir + "mro_sclkscet_00021_65536.tsc")

# Load initial MRO states and store it
MROpos_SSB_ini, MROvel_SSB_ini = spice.get_state("MRO", DYN_TIME_out["time_UTC"])
MRO_out = {
    "MROpos_SSB" : MROpos_SSB_ini,
    "MROvel_SSB" : MROvel_SSB_ini,
}

DYN_SUN_out   = DYN_SUN.initialize(DYN_TIME_out)
DYN_EARTH_out = DYN_EARTH.initialize(DYN_TIME_out)
DYN_MARS_out  = DYN_MARS.initialize(DYN_TIME_out)
DYN_TRA_out   = DYN_TRA.initialize(DYN_EARTH_out, DYN_MARS_out, DYN_SUN_out)
DYN_GRV_out   = DYN_GRV.initialize(DYN_TRA_out)

DYN_out = {
    "DYN_TIME"  : DYN_TIME_out,
    "DYN_SUN"   : DYN_SUN_out,
    "DYN_EARTH" : DYN_EARTH_out,
    "DYN_MARS"  : DYN_MARS_out,
    "DYN_TRA"   : DYN_TRA_out,
    "DYN_GRV"   : DYN_GRV_out,
    "MRO"       : MRO_out,
}

modules = {
    "DYN_TIME"  : DYN_TIME,
    "DYN_SUN"   : DYN_SUN,
    "DYN_EARTH" : DYN_EARTH,
    "DYN_MARS"  : DYN_MARS,
    "DYN_TRA"   : DYN_TRA,
    "DYN_GRV"   : DYN_GRV,
}

# Initialize timeline
timeline = PPC.init_timeline({"DYN" : DYN_out}, n_steps)

# Test loop
for step in range(n_steps):
    # Update algebraic modules first
    DYN_out = DYN_TIME.outputs(sim_time,  DYN_out)
    DYN_out = DYN_SUN.outputs(sim_time,   DYN_out)
    DYN_out = DYN_EARTH.outputs(sim_time, DYN_out)
    DYN_out = DYN_MARS.outputs(sim_time,  DYN_out)
    DYN_out = DYN_TRA.outputs(sim_time,   DYN_out)
    DYN_out = DYN_GRV.outputs(sim_time,   DYN_out)

    # Update MRO state
    MROpos_SSB, MROvel_SSB = spice.get_state("MRO", DYN_out["DYN_TIME"]["time_UTC"])
    DYN_out["MRO"]["MROpos_SSB"] = MROpos_SSB
    DYN_out["MRO"]["MROvel_SSB"] = MROvel_SSB

    # Integrate all dynamic states together
    DYN_out = integrator.rk4_step(sim_time, sim_dt, DYN_out, modules)

    # Save results into the timeline
    output = {"DYN" : DYN_out}
    PPC.update_timeline(timeline, output, step)

    # Update time
    sim_time += sim_dt

# Clear Kernels from memory
spice.clear_kernels()

###########
# RESULTS #
###########

# Plot trajectory
fig, ax = PPC.plot([0], [0], [0], style='+', label="SSB", xlabel="X SSB [km]", ylabel="Y SSB [km]", zlabel="Z SSB [km]", title="Simulation Trajectories")
fig, ax = PPC.plot(timeline["DYN"]["DYN_SUN"]["SUNpos_SSB"][:,0],     timeline["DYN"]["DYN_SUN"]["SUNpos_SSB"][:,1],     timeline["DYN"]["DYN_SUN"]["SUNpos_SSB"][:,2],     label="Sun",   fig=fig, ax=ax)
fig, ax = PPC.plot(timeline["DYN"]["DYN_EARTH"]["EARTHpos_SSB"][:,0], timeline["DYN"]["DYN_EARTH"]["EARTHpos_SSB"][:,1], timeline["DYN"]["DYN_EARTH"]["EARTHpos_SSB"][:,2], label="Earth", fig=fig, ax=ax)
fig, ax = PPC.plot(timeline["DYN"]["DYN_EARTH"]["MOONpos_SSB"][:,0],  timeline["DYN"]["DYN_EARTH"]["MOONpos_SSB"][:,1],  timeline["DYN"]["DYN_EARTH"]["MOONpos_SSB"][:,2],  label="Moon",  fig=fig, ax=ax)
fig, ax = PPC.plot(timeline["DYN"]["DYN_MARS"]["MARSpos_SSB"][:,0],   timeline["DYN"]["DYN_MARS"]["MARSpos_SSB"][:,1],   timeline["DYN"]["DYN_MARS"]["MARSpos_SSB"][:,2],   label="Mars",  fig=fig, ax=ax)
fig, ax = PPC.plot(timeline["DYN"]["DYN_MARS"]["DEIMOSpos_SSB"][:,0], timeline["DYN"]["DYN_MARS"]["DEIMOSpos_SSB"][:,1], timeline["DYN"]["DYN_MARS"]["DEIMOSpos_SSB"][:,2], label="Deimos",fig=fig, ax=ax)
fig, ax = PPC.plot(timeline["DYN"]["DYN_MARS"]["PHOBOSpos_SSB"][:,0], timeline["DYN"]["DYN_MARS"]["PHOBOSpos_SSB"][:,1], timeline["DYN"]["DYN_MARS"]["PHOBOSpos_SSB"][:,2], label="Phobos",fig=fig, ax=ax)
fig, ax = PPC.plot(timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,0],      timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,1],      timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,2],      label="SC",    fig=fig, ax=ax)
fig, ax = PPC.plot(timeline["DYN"]["MRO"]["MROpos_SSB"][:,0],         timeline["DYN"]["MRO"]["MROpos_SSB"][:,1],         timeline["DYN"]["MRO"]["MROpos_SSB"][:,2],         label="MRO",   fig=fig, ax=ax)

# Compare MRO with SC
fig, ax_x = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["MRO"]["MROpos_SSB"][:,0], label="MROx_SSB", ylabel="[km]", subplot=[4,1,1])
PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,0], label="SCx_SSB", fig=fig, ax=ax_x)
fig, ax_y = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["MRO"]["MROpos_SSB"][:,1], label="MROy_SSB", ylabel="[km]", fig=fig, subplot=[4,1,2])
PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,1], label="SCy_SSB", fig=fig, ax=ax_y)
fig, ax_z = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["MRO"]["MROpos_SSB"][:,2], label="MROz_SSB", xlabel="time_SIM [s]", ylabel="[km]", fig=fig, subplot=[4,1,3])
PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,2], label="SCz_SSB", fig=fig, ax=ax_z)

error_x =  timeline["DYN"]["MRO"]["MROpos_SSB"][:,0] - timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,0]
error_y =  timeline["DYN"]["MRO"]["MROpos_SSB"][:,1] - timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,1]
error_z =  timeline["DYN"]["MRO"]["MROpos_SSB"][:,2] - timeline["DYN"]["DYN_TRA"]["SCpos_SSB"][:,2]
total_error = (error_x**2 + error_y**2 + error_z**2)**0.5
fig, ax_error = PPC.plot(timeline["DYN"]["DYN_TIME"]["time_SIM"], total_error, label="SSB-MRO error", fig=fig, subplot=[4,1,4])

PPC.show_plot()