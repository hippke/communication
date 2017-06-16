import math
from astropy.constants import c, G, M_sun, R_sun, au, L_sun, pc
from astropy import units as u
from math import pi, sqrt, exp, log, log2
from PyCom import holevo_perfect, holevo_thermal, photons_received, b_of_z, \
    z_of_b, noise_photons, sgl_psf, integrated_flux, classical_aperture, d_critical


"""Calculates all numbers used in part 2 of the series"""


def angular_resolution(d, wavelength, Q=1.22):
    """Returns the angular resolution of a telescope [arcsec]
    d: circular aperture diameter [m]
    wavelength \lambda [m]
    Q: Quality < 1.22, diffraction-limit: Q=1.22"""
    return 3600 * 180 / pi * Q * wavelength / d


def apparent_size_object(distance, size_object):
    """Returns the apparent diameter of an object in [arcsec]"""
    return 3600 * 180 / pi * 2 * size_object / u.meter / distance


### Section 2.1 Lens geometry

# Schwarzschild radius of the sun
r_g_sun = ((2 * G * M_sun) / (c**2)) / u.meter  # [m] Schwarzschild radius of the sun
print('Schwarzschild radius of the sun [m]', '{:.2f}'.format(r_g_sun))

# Figure caption: How many R_sun are 546 au?
print('546 au are {:.1e} R_sun'.format(546 * au / R_sun))

# z(b) and b(z)
print('z(b=1.05) [au]', '{:.2f}'.format(z_of_b(b=1.05) / au * u.meter))
print('b(z=1000 au)', '{:.2f}'.format(b_of_z(z=1000 * au / u.meter)))

# Apparent width of Einstein ring
d = 1  # [m] aperture of SGL telescope
se = apparent_size_object(distance=600 * au / u.meter, size_object=d / 2)
print('Apparent size of Einstein ring [arcsec]', '{:.1e}'.format(se))

# Angular resolution of telescope
wavelength = 1000 * 10**-9  # [m] 1000 nm = 1um = 10**-6m
ar = angular_resolution(d=d, wavelength=wavelength)
print('Angular resolution', '{:.2f}'.format(ar))

# Ratio of Einstein ring and angular resolution (how much is it unresolved?)
print('Einstein ring unresolved by', '{:.2e}'.format(se / ar))


### Section 2.2 Frequency cutoff
frequency = 410  # [GHz]
critical_distance = d_critical(frequency=frequency) / au * u.meter
print('Minimum focusing distance for frequency', frequency, '[GHz]: {:.1f} [au]'.format(critical_distance))


### Section 2.3 Light collection power

# Gain at z=600
eq8_at_z_600 = integrated_flux(wavelength=wavelength, z=600 * au / u.meter, d0=d)
print('Gain at z=600 {:.2e}'.format(eq8_at_z_600))

# Move it to z=2200 to compare influence of z
eq8_at_z_2200 = integrated_flux(wavelength=wavelength, z=2200 * au / u.meter, d0=d)
print('Gain at z=2200 {:.2e}'.format(eq8_at_z_2200))
print('Ratio z=600 / z=2200 {:.2f}'.format(eq8_at_z_2200 / eq8_at_z_600))

# Equivalent classical aperture
ap_600 = classical_aperture(wavelength=wavelength, z=600 * au / u.meter, d0=d)
print('Classical aperture at z=600 {:.2f} [km]'.format(ap_600 / 1000))


### Section 3.1 Photon flux
flux = photons_received(D_r=1, D_t=1, P=1, wavelength=wavelength, R=1.3 * pc / u.meter, Q_R=1.22)
print('Signal photons received (d_r=d_t=1m, P=1W, Alpha Cen) {:.2e} [m^-2 s^-1]'.format(flux))


### Section 4.6 Noise calculations
noise = noise_photons(
    wavelength=wavelength,
    z=600 * au / u.meter,
    d0=d,
    fractional_bandwidth=1/2700.)
print('Noise photons received {:.2e} [s^-1]'.format(noise))


### Section 5.1 Capacity

receiver_efficiency = 0.5
F = 1  # number of photons
N = 0.13  # Number of noise photons
modes = 10**5  # number of modes
M = F / modes  # Number of photons per mode (the signal)
N_th = 0.13  #  N / modes  # average number of noise photons per mode (the noise)
capacity = holevo_thermal(M=M, N_th=N_th, eta=receiver_efficiency)
print('Capacity {:.2f} [bits per photon]'.format(capacity))


### Section 5.2 Data rate
# Signal
wavelength = 1000 * 10**-9  # 1000 nm = 1um = 10**-6m
result1 = photons_received(D_r=1, D_t=1, P=1, wavelength=wavelength, R=1.3 * pc / u.meter, Q_R=1.22)
print('A pair of classical 1m telescopes at a distance of 1.3 pc transmitting at\
 lambda=1mu can exchange', '{:.2e}'.format(float(result1)), 'photons per Watt per second, neglecting losses.')

result2 = classical_aperture(wavelength=0.000001, z=600 * au / u.meter, d0=1)
print('The same 1m receiving telescope in the SGL with its equivalent aperture of',\
    '{:.2f}'.format(float(result2 / 1000)), 'km')

result3 = photons_received(D_r=result2, D_t=1, P=1, wavelength=wavelength, R=1.3 * pc / u.meter, Q_R=1.22)
print('collects', '{:.2e}'.format(float(result3)), 'photons per second per Watt, \
a gain of', '{:.2e}'.format(float(result3 / result1)))

result4 = classical_aperture(wavelength=0.000001, z=2200 * au / u.meter, d0=1)
result5 = photons_received(D_r=result4, D_t=1, P=1, wavelength=wavelength, R=1.3 * pc / u.meter, Q_R=1.22)
print('At z=600 the equivalent D_classical is', '{:.2f}'.format(float(result4 / 1000)), 'km')
print('For a flux of ', '{:.2e}'.format(float(result5)), 'photons per Watt per second, neglecting losses.')

wavelength = 1000 * 10**-9  # 1000 nm = 1um = 10**-6m
diameter = 1
fractional_bandwidth = 1/2700
power = 1  # Watt
modes = 10**10
z = 600  # au

# Number of noise photons
noise = noise_photons(
    wavelength=wavelength,
    z=z*au/u.meter,
    d0=diameter,
    fractional_bandwidth=fractional_bandwidth)

# Classical aperture equivalent
aperture = classical_aperture(
    wavelength=wavelength,
    z=z*au/u.meter,
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
print('Data rate d_SGL, 1 meter', data_rate)

signal = photons_received(
    D_r=49000,
    D_t=1,
    P=power,
    wavelength=wavelength,
    R=1.3 * pc / u.meter,
    Q_R=1.22)
M = signal / modes  # Number of photons per mode (the signal)
N_th = 0.1 / modes  # average number of noise photons per mode (the noise)
capacity = holevo_thermal(M=M, N_th=N_th, eta=receiver_efficiency)
data_rate = (capacity * signal) / 10**6  # Mbits/s
print('Data rate d_classical, 45 km', data_rate)

a1 = 49  # [km] aperture
a2 = 53.57  # [km] aperture

a1_area = pi * a1**2 / 4
a2_area = pi * a2**2 / 4
realized_aperture_fraction = a1_area / a2_area
print('Realized aperture fraction', realized_aperture_fraction)