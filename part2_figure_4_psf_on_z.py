import numpy
from astropy import units as u
from astropy.constants import c, G, M_sun, R_sun, au
from matplotlib import pyplot as plt
from PyCom import integrated_flux, b_of_z


wavelength = 10**-6  # [m]
path = ''
list_b0 = []
list_flux = []
new_d0 = 1
base_flux = integrated_flux(wavelength=wavelength, z=(550 * au) / u.meter, d0=new_d0)
for i in range(546, 2200):
    z = (i * au) / u.meter
    new_flux = integrated_flux(wavelength=wavelength, z=z, d0=new_d0)
    b = b_of_z(z)
    list_b0.append(b)
    list_flux.append(new_flux / base_flux)

plt.rc('font',  family='serif', serif='Computer Modern Roman')
plt.rc('text', usetex=True)
size = 4.5
aspect_ratio = 1.5
plt.figure(figsize=(size, size / aspect_ratio))
ax = plt.gca()
ax.get_yaxis().set_tick_params(direction='out')
ax.get_xaxis().set_tick_params(direction='out')
ax.get_yaxis().set_tick_params(which='both', direction='out')
ax.get_xaxis().set_tick_params(which='both', direction='out')
ax.set_xlabel(r'$R$ $(R_{\odot})$')
ax.set_ylabel(r'Relative flux')
ax.set_xlim(1, 2)
ax.set_ylim(1, 2)
plt.plot(list_b0, list_flux, color='black')
ay2 = ax.twiny()
ay2.get_xaxis().get_major_formatter().set_useOffset(False)
ay2.get_xaxis().get_major_formatter().set_scientific(False)
ay2.get_xaxis().set_tick_params(direction='out')
ay2.get_xaxis().set_tick_params(direction='out')
ay2.get_xaxis().set_tick_params(which='both', direction='out')
ay2.get_xaxis().set_tick_params(which='both', direction='out')
ay2.set_xlim([546, 2184])
ay2.set_xlabel(r'Heliocentric distance $z$ (au)')
plt.savefig(path + 'figure_z.pdf', bbox_inches='tight')
