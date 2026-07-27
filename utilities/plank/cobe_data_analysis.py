import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import healpy as hp
from astropy.io import fits
from scipy.spatial import cKDTree

# COBE (not in HEALPix format)
# wget https://lambda.gsfc.nasa.gov/data/cobe/dmr/asds4/cosmic_emission_maps/DMR_DCMB_GALACTIC_4YR.FITS
# wget https://lambda.gsfc.nasa.gov/data/cobe/dmr/asds4/cosmic_emission_maps/DMR_DSMB_GALACTIC_4YR.FITS
filename = "DMR_DCMB_GALACTIC_4YR.FITS"; scale = 1e-3 # [mK]

# 1. Load data and coordinates from FITS
with fits.open(filename) as hdul:
    table = hdul['DMR_SKYMAP'].data
    signal = table['SIGNAL']
    l_deg = table['GALON']
    b_deg = table['GALAT']

# 2. Convert COBE Galactic coordinates to 3D Cartesian vectors
l_rad = np.radians(l_deg)
b_rad = np.radians(b_deg)
x_cobe = np.cos(b_rad) * np.cos(l_rad)
y_cobe = np.cos(b_rad) * np.sin(l_rad)
z_cobe = np.sin(b_rad)
cobe_vecs = np.column_stack([x_cobe, y_cobe, z_cobe])

# 3. Build 3D spatial search tree for COBE pixels
tree = cKDTree(cobe_vecs)

# 4. Generate 3D vectors for ALL HEALPix pixel centers (NSIDE=32)
nside = 64
npix = hp.nside2npix(nside)
hpx_vecs = np.array(hp.pix2vec(nside, np.arange(npix))).T

# 5. Query nearest COBE pixel for EVERY target HEALPix pixel (Eliminates empty pixel holes)
_, nearest_indices = tree.query(hpx_vecs)
cmb_map = signal[nearest_indices] * scale

# 6. Smooth with COBE's 7-degree FWHM instrument beam
cmb_map_smoothed = hp.smoothing(cmb_map, fwhm=np.radians(7.0))

cmb_cmap = LinearSegmentedColormap.from_list(
    "cmb",
    [
        (0.00,(0.00, 0.40, 1.00)), # blue
        # (0.20, 0.80, 0.00), # green
        (0.5,(1,1,1)), # white
        (1.00,(0.80, 0.00, 0.00)), # red
    ],
    N=512,
)

# 7. Plot continuous sky map
hp.mollview(
    cmb_map_smoothed,
    title="COBE DMR 4-Year Galactic Emission Map",
    unit="K",
    cmap=cmb_cmap,
    coord="G",
    min=-300e-6,
    max=300e-6
)
hp.graticule()

plt.show()

###############################
# EXPORT VECTOR HAMMER MAP    #
###############################

plt.figure(figsize=(12, 6))

hp.projview(
    cmb_map_smoothed,
    projection_type="hammer",
    coord=["G"],            # Galactic coordinates
    cmap=cmb_cmap,
    graticule=False,
    cbar=False,
    title="",
    xlabel="",
    ylabel="",
    flip="astro",           # Astronomical convention (l increases to the left)
    min=-300e-6,
    max=300e-6
)

# Remove any remaining axes decorations
ax = plt.gca()
ax.set_axis_off()

plt.savefig(
    "cmb_cobe.pdf",
    format="pdf",
    transparent=True,
    bbox_inches="tight",
    pad_inches=0,
)

plt.close()