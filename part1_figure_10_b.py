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
D_r = 39  # m
D_t = 1  # m
P = 1000  # Watt
eta = 0.5
N = 10**-5  # modes per photon ^-1

gridsize = 500
image = numpy.zeros(shape=(gridsize,gridsize))
distances = numpy.logspace(start=0, stop=4, num=gridsize)
wavelengths = numpy.logspace(start=2, stop=4, num=gridsize)

dist_id = 0
wave_id = 0

for distance in distances:
    dist_id = dist_id + 1
    wave_id = 0
    for wavelength in wavelengths:
        wave_id = wave_id + 1
        p = photons_received(D_r=D_r, D_t=D_t, P=P, wavelength=wavelength, R=distance * parsec, Q_R=Q(wavelength))
        e = get_extinction(distance=distance, wavelength=wavelength)
        tot = p * e #* t * bits_per_photon
        image[dist_id-1,wave_id-1] = tot  # log(1/tot)

    print(distance)


for row in range(gridsize):
    image[row,:] = image[row,:] / numpy.amax(image[row,:])

plt.rc('font',  family='serif', serif='Computer Modern Roman')
plt.rc('text', usetex=True)
size = 7.5
aspect_ratio = 1.5
plt.figure(figsize=(size, size / aspect_ratio))
ax = plt.gca()
i = ax.imshow(image, interpolation='nearest', origin='lower', cmap=plt.get_cmap('inferno'))
#colorbar_ax = fig.add_axes([0.7, 0.1, 0.05, 0.8])
#fig.colorbar(i, cax=colorbar_ax)
ax.set_xticks((0, int(gridsize/2), gridsize-1))
ax.set_xticklabels(('0.1', '1', '10'))
ax.set_yticks((0, int(1*gridsize/4), int(2*gridsize/4), int(3*gridsize/4), gridsize-1))
ax.set_yticklabels(('1', '10', '100', '1000', '10000'))
ax.get_yaxis().set_tick_params(direction='out')
ax.get_xaxis().set_tick_params(direction='out')
ax.get_yaxis().set_tick_params(which='both', direction='out')
ax.get_xaxis().set_tick_params(which='both', direction='out')
ax.set_xlabel(r'Wavelength $\lambda$ ($\mu$m)')
ax.set_ylabel(r'Distance (pc)')
plt.savefig(path + 'figure_grid_photons_space.pdf', bbox_inches='tight')
#plt.show()
