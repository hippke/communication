from PyCom import holevo_perfect, holevo_thermal, photons_received, Q
import numpy
from math import log, log2, pi
from matplotlib import pyplot as plt
import matplotlib.patches as patches

path = ''


def sky_background(wavelength):
    """Returns the sky background in phot/s/nm/arcsec^2/m^2 for wavelength
    Wavelength must be in nm"""
    idx = numpy.argmin(numpy.abs(data_sky_background['wavelength'] - wavelength))
    return data_sky_background['flux'][idx]


def transparency(wavelength):
    """Returns the sky transparency in [0, 1] for wavelength (in nm)"""

    idx = numpy.argmin(numpy.abs(
        data_transparency['wavelength'] * 1000 - wavelength))
    return data_transparency['fraction'][idx]


def get_extinction(distance, wavelength):
    scale_factor = 0.2665
    A = draine['C_ext'] / scale_factor * 10**21  # normalized to 1000 pc
    curve = 1 / 2.512**(A/(1000/distance))
    curve = curve / lyman_extinct['factor']
    idx = numpy.argmin(numpy.abs(draine['wavelength'] - wavelength / 1000))
    return curve[idx]


draine = numpy.genfromtxt(
    path + 'data_kext_albedo_WD_MW_3.1_60_D03.all',
    skip_header=80,
    usecols=(0,1,2,3,4,5),
    autostrip=True,
    dtype=[
        ('wavelength', 'f8'),
        ('albedo', 'f8'),
        ('cos', 'f8'),
        ('C_ext', 'f8'),
        ('K_abs', 'f8'),
        ('cos2', 'f8')])

lyman_extinct = numpy.genfromtxt(
    path + 'data_lyman_extinct.csv',
    dtype=[('factor', 'f8')])

data_sky_background = numpy.genfromtxt(path + 'data_sky_background_collated.txt',
    dtype=[
        ('wavelength', 'f8'),
        ('flux', 'f8')])

data_transparency = numpy.genfromtxt(path + 'data_transparency.txt',
    dtype=[
        ('wavelength', 'f8'),
        ('fraction', 'f8')])


parsec = 3.086e+16  # meter

wavelength = 429.6  # nm
distance = 1.3  # pc
D_r = 10  # m
D_t = 1  # m
P = 1000  # Watt
eta = 0.5
N = 10**-5  # modes per photon ^-1

list_waves = []
list_photons = []

for step in range(1, 4000):
    wavelength = step
    p = photons_received(D_r=D_r, D_t=D_t, P=P, wavelength=wavelength, R=distance * parsec, Q_R=Q(wavelength))
    e = get_extinction(distance=distance, wavelength=wavelength)
    receiver_aperture = pi * (D_t / 2)**2
    Noise = 0.1 * receiver_aperture
    N_th = N * Noise  # Noise photons per mode
    bits_per_photon = holevo_thermal(N=N, N_th=N_th, eta=eta)
    tot = p * e * bits_per_photon
    print(wavelength, Noise, bits_per_photon, tot, p, e)
    list_waves.append(wavelength/1000)
    list_photons.append(tot)


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
ax.set_xlabel(r'Wavelength $\lambda$ ($\mu$m)')
ax.set_ylabel(r'Data rate (bits/s)')
ax.set_xlim(0.1, 4)
ax.set_ylim(0, 900)
plt.plot(list_waves, list_photons, color='black', linewidth=0.5)

cm = plt.cm.get_cmap('jet')
for i in range(1000):
    ax.add_patch(patches.Rectangle((0.35 + i/1000/3, 0), 0.01, 60,
        facecolor=cm(i/1000), edgecolor='none'))

ax.text(0.72, 10, r'IR')
ax.text(3.9, 810, r'$D_t=1$\,m (space), $D_r=10$\,m (space)', horizontalalignment='right')
ax.text(3.9, 725, r'$R=1.3$\,pc, $P=1000$\,W', horizontalalignment='right')
ax.text(3.9, 640, r'$\eta=0.5$, $M=10^5$', horizontalalignment='right')
plt.savefig(path + 'proxima_space.pdf', bbox_inches='tight')
#plt.show()
