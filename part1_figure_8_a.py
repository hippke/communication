from PyCom import holevo_perfect, holevo_thermal, photons_received, Q
import numpy
from math import log, log2, pi
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from astropy.constants import c, G, M_sun, R_sun, au, L_sun, pc
from astropy import units as u


path = ''


def sky_background(wavelength):
    """Returns the sky background in phot/s/nm/arcsec^2/m^2 for wavelength
    Wavelength must be in [m]"""
    wavelength = wavelength / 10**-9
    idx = numpy.argmin(numpy.abs(data_sky_background['wavelength'] - wavelength))
    return data_sky_background['flux'][idx]


def transparency(wavelength):
    """Returns the sky transparency in [0, 1] for wavelength in [m]"""
    wavelength = wavelength / 10**-9
    idx = numpy.argmin(numpy.abs(
        data_transparency['wavelength'] * 1000 - wavelength))
    return data_transparency['fraction'][idx]


def get_extinction(distance, wavelength):
    """distance in [pc]"""
    scale_factor = 0.2665
    wavelength = wavelength / 10**-9
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


distance = 1.3  # pc
D_r = 39  # m
D_t = 1  # m
P = 1000  # Watt
eta = 0.5
modes = 10**5  # modes per photon ^-1

list_waves = []
list_photons = []

for step in range(2900, 40000):
    wavelength = step * 10**-9 / 10 # [m] to [nm], finer sampling * 10
    p = photons_received(D_r=D_r, D_t=D_t, P=P, wavelength=wavelength, R=distance * pc / u.meter, Q_R=Q(wavelength))
    e = get_extinction(distance=distance, wavelength=wavelength)
    t = transparency(wavelength=wavelength)
    F = p * e * t
    receiver_aperture = pi * (D_r / 2)**2
    N = sky_background(wavelength=wavelength) * receiver_aperture
    N_th = N / modes  # average number of noise photons per mode (the noise)
    M = F / modes  # Number of photons per mode (the signal)
    bits_per_photon = holevo_thermal(M=M, N_th=N_th, eta=eta)
    tot = F * bits_per_photon
    print(wavelength * 10**6, F, N, bits_per_photon, p, e, t, tot)
    list_waves.append(wavelength * 10**6)
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
ax.set_ylim(0, 2600)
plt.plot(list_waves, list_photons, color='black', linewidth=0.5)

cm = plt.cm.get_cmap('jet')
for i in range(1000):
    ax.add_patch(patches.Rectangle((0.35 + i/1000/3, 0), 0.01, 200,
        facecolor=cm(i/1000), edgecolor='none'))
ax.text(0.72, 30, r'IR')
ax.text(3.9, 2300, r'$D_{\rm t}=1$\,m (space), $D_{\rm r}=39$\,m (Earth)', horizontalalignment='right')
ax.text(3.9, 2000, r'$R=1.3$\,pc, $P=1000$\,W', horizontalalignment='right')
ax.text(3.9, 1700, r'$\eta=0.5$, $M=10^5$', horizontalalignment='right')
plt.savefig(path + 'proxima_earth.pdf', bbox_inches='tight')
