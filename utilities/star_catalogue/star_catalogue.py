# Loads Hipparcos star catalogue using skyfield library

from skyfield.api import load, Star
from skyfield.data import hipparcos
import matplotlib.pyplot as plt
import numpy as np

# Load entire star catalogue
with load.open("hip_main.dat") as f:
    star_df = hipparcos.load_dataframe(f)

# Filter stars by magnitude
star_df = star_df[star_df["magnitude"] <= 3]

stars = Star.from_dataframe(star_df)
pos = stars._position_au
norm = np.linalg.norm(pos, axis=0)
x, y, z = pos/norm

# 2D plot
plt.scatter(star_df["ra_degrees"], star_df["dec_degrees"],marker='.')
plt.xlabel('RA [deg]')
plt.ylabel('Decl [deg]')
plt.grid()

# 3D plot
fig = plt.figure()
ax  = fig.add_subplot(1,1,1, projection='3d')
ax.plot(x, y, z, '.')
ax.set_xlim([-1, 1])
ax.set_ylim([-1, 1])
ax.set_zlim([-1, 1])
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.set_aspect("equal")

# Draw celestial sphere outline
u = np.linspace(0, 2*np.pi, 100)
v = np.linspace(-np.pi/2, np.pi/2, 100)
xs = np.cos(v)[:, None] * np.cos(u)
ys = np.cos(v)[:, None] * np.sin(u)
zs = np.sin(v)[:, None]
ax.plot_wireframe(xs, ys, zs, color='lightgray', alpha=0.2)

plt.show()
