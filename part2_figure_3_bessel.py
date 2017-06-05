import numpy
from math import pi, sqrt, exp
from scipy.special import j0, j1  # Bessel function of the first kind of order 0
from astropy.constants import c, G, M_sun, R_sun, au
from astropy import units as u
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from PyCom import sgl_psf, classical_aperture


wavelength = 10**-6  # [m]
z = (600 * au) / u.meter  # [m] heliocentric distance from the sun

# Figure 3a: SGL PSF
array_size = 10000
rho_array = numpy.linspace(0, 10, array_size)
flux_array = numpy.linspace(0, 10, array_size)
for value in range(len(rho_array)):
    flux_array[value] = sgl_psf(rho=value/1000, wavelength=wavelength, z=z)

plt.rc('font',  family='serif', serif='Computer Modern Roman')
plt.rc('text', usetex=True)
size = 4.5
aspect_ratio = 1.5
plt.figure(figsize=(size, size / aspect_ratio))
ax = plt.gca()
ax.get_yaxis().get_major_formatter().set_useOffset(False)
ax.get_yaxis().get_major_formatter().set_scientific(False)
ax.get_yaxis().set_tick_params(direction='out')
ax.get_xaxis().set_tick_params(direction='out')
ax.get_yaxis().set_tick_params(which='both', direction='out')
ax.get_xaxis().set_tick_params(which='both', direction='out')
ax.set_ylim([0, 1.05])
ax.set_xlim([-0.3, 0.3])
ax.set_xlabel(r'Distance from axis (m)')
ax.set_ylabel(r'Relative flux')
plt.plot(rho_array, flux_array / flux_array[0], color='black')
plt.plot(-rho_array, flux_array / flux_array[0], color='black')
plt.savefig('bessel_a.pdf', bbox_inches='tight')

# Figure 3b: Classical versus SGL telescope size
list_d0 = []
list_flux = []
list_size = []
for i in range(1, 10000):
    new_d0 = i / 100
    new_classical_aperture = classical_aperture(wavelength=wavelength, z=z, d0=new_d0) / 1000
    list_d0.append(new_d0)
    list_size.append(new_classical_aperture)

plt.figure(figsize=(size, size / aspect_ratio))
ax = plt.gca()
ax.get_yaxis().get_major_formatter().set_useOffset(False)
ax.get_yaxis().get_major_formatter().set_scientific(False)
ax.get_yaxis().set_tick_params(direction='out')
ax.get_xaxis().set_tick_params(direction='out')
ax.get_yaxis().set_tick_params(which='both', direction='out')
ax.get_xaxis().set_tick_params(which='both', direction='out')
ax.set_xlabel(r'SGL telescope diameter (m)')
ax.set_ylabel(r'Classical telescope diameter (km)')
plt.plot(list_d0, list_size, color='black')
ax.set_xlim([0, 1])
ax.set_ylim([0, 60])
plt.savefig('telescope_size_a.pdf', bbox_inches='tight')


# Figure 3c: ZOOM for classical versus SGL telescope size
ax.set_xlim([0, 100])
ax.set_ylim([0, 600])
ax.add_patch(patches.Rectangle((0, 0), 1, 60,))
plt.savefig('telescope_size_b.pdf', bbox_inches='tight')
