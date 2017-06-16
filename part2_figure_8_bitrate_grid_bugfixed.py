import numpy
from astropy.constants import c, G, M_sun, R_sun, au, L_sun, pc
from astropy import units as u
from math import pi, sqrt, exp, log, log2
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from PyCom import classical_aperture, noise_photons, holevo_thermal, photons_received


# Astrophysical parameters
modes = 10**10  # number of modes
receiver_efficiency = 0.5
wavelength=0.000001
fractional_bandwidth = 1/2700
power = 1  # Watt

# Figure parameters
gridsize = 250
image = numpy.zeros(shape=(gridsize,gridsize))
distances = numpy.linspace(start=547, stop=2200, num=gridsize)
diameters = numpy.linspace(start=0.01, stop=12, num=gridsize)
dist_id = 0
diam_id = 0

for distance in distances:
    dist_id = dist_id + 1
    diam_id = 0
    for diameter in diameters:
        diam_id = diam_id + 1

        # Number of noise photons
        noise = noise_photons(
            wavelength=wavelength,
            z=distance*au/u.meter,
            d0=diameter,
            fractional_bandwidth=fractional_bandwidth)

        # Classucal aperture equivalent
        aperture = classical_aperture(
            wavelength=wavelength,
            z=distance*au/u.meter,
            d0=diameter)

        # number of signal photons
        signal = photons_received(
            D_r=aperture,
            D_t=1,
            P=power,
            wavelength=wavelength,
            R=1.3 * pc / u.meter,
            Q_R=1.22)

        M = signal / modes  # Number of photons per mode (the signal)
        N_th = noise / modes  # average number of noise photons per mode (the noise)
        capacity = holevo_thermal(M=M, N_th=N_th, eta=receiver_efficiency)

        SNR = signal / noise
        data_rate = (capacity * signal) / 10**6  # Mbits/s
        image[dist_id-1,diam_id-1] = data_rate
        # print(diameter, distance, SNR, capacity, data_rate)
    print(distance)

print(numpy.amax(image), numpy.amin(image))
plt.rc('font',  family='serif', serif='Computer Modern Roman')
plt.rc('text', usetex=True)
size = 5.5
aspect_ratio = 1.5
fig, ax = plt.subplots(figsize=(size, size / aspect_ratio))
i = ax.imshow(image, interpolation='nearest', origin='lower',
    cmap=plt.get_cmap('inferno'), vmin=0, vmax=93)  # inferno jet magma, norm=LogNorm(vmin=1, vmax=100)
colorbar_ax = fig.add_axes()  # [0.7, 0.1, 0.05, 0.8]
fig.colorbar(i, cax=colorbar_ax)
ticks = (0, int(1*gridsize/4), int(2*gridsize/4), int(3*gridsize/4), gridsize-1)
ax.set_xticks(ticks)
ax.set_xticklabels(('0', '3', '6', '9', '12'))
ax.set_yticks(ticks)
ax.set_yticklabels(('547', '960', '1374', '1787', '2200'))
ax.get_yaxis().set_tick_params(direction='out')
ax.get_xaxis().set_tick_params(direction='out')
ax.get_yaxis().set_tick_params(which='both', direction='out')
ax.get_xaxis().set_tick_params(which='both', direction='out')
ax.set_xlabel(r'SGL aperture $d$ (m)')
ax.set_ylabel(r'Heliocentric distance $z$ (au)')
plt.savefig('figure_sgl_bitrate_grid.pdf', bbox_inches='tight')
#plt.show()
