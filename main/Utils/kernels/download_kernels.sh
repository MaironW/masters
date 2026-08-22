# Kernel files are too large to store using git
# Run this script to get the expected kernels to run the simulations

wget https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls
wget https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp
wget https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc
wget https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/a_old_versions/mar097.bsp
wget https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_000101_260114_251018.bpc

# MRO kernel for tests
wget https://naif.jpl.nasa.gov/pub/naif/pds/data/mro-m-spice-6-v1.0/mrosp_1000/data/spk/mro_cruise.bsp

# MSL kernal for tests
wget https://naif.jpl.nasa.gov/pub/naif/MSL/kernels/spk/msl_cruise.bsp

# Not a Kernel, but star catalogue
wget https://cdsarc.cds.unistra.fr/ftp/cats/I/239/hip_main.dat

# CMBR WMAP data (because it is lighter than the Plank one)
wget https://lambda.gsfc.nasa.gov/data/map/dr5/dfp/ilc/wmap_ilc_9yr_v5.fits