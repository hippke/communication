import math
from astropy.constants import c, G, M_sun, R_sun, au, L_sun, pc
from astropy import units as u
from math import pi, sqrt, exp, log, log2
from PyCom import holevo_perfect, holevo_thermal, photons_received, b_of_z, \
    z_of_b, noise_photons, sgl_psf, integrated_flux, classical_aperture


"""Test Holevo thermal"""
receiver_efficiency = 0.5

# Case 1
F = 100  # number of photons
N = 100  # Number of noise photons
modes = 10**5  # number of modes
M = F / modes  # Number of photons per mode (the signal)
N_th = N / modes  # average number of noise photons per mode (the noise)
result = holevo_thermal(M=M, N_th=N_th, eta=receiver_efficiency)
expected_result = 2.6023902315676968
if math.isclose(result, expected_result, rel_tol=1e-6):
    print('PASS holevo_thermal test case 1, result', result)
else:
    print('FAILED holevo_thermal test case 1, result', result)

# Case 2
F = 50  # number of photons
N = 10  # Number of noise photons
M = F / modes  # Number of photons per mode (the signal)
N_th = N / modes  # average number of noise photons per mode (the noise)
result = holevo_thermal(M=M, N_th=N_th, eta=receiver_efficiency)
expected_result = 3.157176216339997
if math.isclose(result, expected_result, rel_tol=1e-6):
    print('PASS holevo_thermal test case 2, result', result)
else:
    print('FAILED holevo_thermal test case 2, result', result)


"""Test Holevo perfect"""
F = 50  # number of photons
M = F / modes  # Number of photons per mode (the signal)
result = holevo_perfect(M=M, eta=receiver_efficiency)
expected_result = 3.352164911851362
if math.isclose(result, expected_result, rel_tol=1e-6):
    print('PASS holevo_perfect, result', result)
else:
    print('FAILED holevo_perfect, result', result)


"""Test z and b"""
result = b_of_z(z=600 * au / u.meter)
expected_result = 1.047036202450992
if math.isclose(result, expected_result, rel_tol=1e-6):
    print('PASS b_of_t test case 1, result', result)
else:
    print('FAILED b_of_t test case 1, result', result)

result = z_of_b(b=2.0) / au * u.meter
expected_result = 2189.2121278750956
if math.isclose(result, expected_result, rel_tol=1e-6):
    print('PASS z_of_b, result', result)
else:
    print('FAILED z_of_b, result', result)


"""Test photons_received"""
wavelength = 1000 * 10**-9  # 1000 nm = 1um = 10**-6m
result = photons_received(D_r=39, D_t=1, P=1, wavelength=wavelength, R=3.09E+16, Q_R=1.22)
expected_result = 1.3469636550722477
if math.isclose(result, expected_result, rel_tol=1e-6):
    print('PASS photons_received, result', result)
else:
    print('FAILED photons_received, result', result)


"""Test noise_photons"""
result = noise_photons(
    wavelength=0.000001,
    z=600 * au / u.meter,
    d0=1,
    fractional_bandwidth=1/2700.)
expected_result = 1905694.5888717826
if math.isclose(result, expected_result, rel_tol=1e-6):
    print('PASS noise_photons, result', result)
else:
    print('FAILED noise_photons, result', result)


"""Test sgl_psf"""
result = sgl_psf(rho=0.0, wavelength=0.000001, z=600 * au / u.meter)
expected_result = 116622066984.28542
if math.isclose(result, expected_result, rel_tol=1e-6):
    print('PASS sgl_psf, result', result)
else:
    print('FAILED sgl_psf, result', result)


"""Test integrated_flux"""
result = integrated_flux(wavelength=0.000001, z=600 * au / u.meter, d0=1)
expected_result = 2869792678.218271
if math.isclose(result, expected_result, rel_tol=1e-6):
    print('PASS integrated_flux, result', result)
else:
    print('FAILED integrated_flux, result', result)
