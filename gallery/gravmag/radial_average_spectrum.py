"""
Radial Average Power Spectrum
-----------------------------

add description

"""
from __future__ import division, print_function
from fatiando import mesher, gridder, utils
from fatiando.datasets import fetch_hawaii_gravity
from fatiando.gravmag import transform
import matplotlib.pyplot as plt
import numpy as np

# Fetch Hawaii data
data = fetch_hawaii_gravity()

# When the data is projected using UTM zone 4, it is no longer regular gridded.
# We must regrid it.
area = (1.483e6, 3.079e6, -60e3, 1.326e6)
shape = (77, 67)
x, y, gravity = gridder.interp(data['x'], data['y'],
                               data['topo-free'],
                               shape, area=area)

# Lets compute the Power Density Spectra (2d arrays)
ky, kx, pds = transform.power_density_spectra(x, y, gravity, shape)

# And then compute the radial average of the PDS
k_radial, pds_radial = transform.radial_average_spectrum(kx, ky, pds)

# Plot results
plt.contourf(y.reshape(shape), x.reshape(shape),
             gravity.reshape(shape), 100)
plt.axes().set_aspect('equal')
plt.colorbar()
plt.xlabel("m")
plt.ylabel("m")
plt.title("Topography Free Gravity Disturbance")
plt.show()

plt.plot(k_radial, np.log(pds_radial), 'o-')
plt.xlabel("log(PDS)")
plt.ylabel("k [1/m]")
plt.show()
