import numpy
from astropy.constants import c, G, M_sun, R_sun, au, L_sun, pc
from astropy import units as u
from math import pi, sqrt, exp, log, log2
from scipy.special import j0, j1  # Bessel function
from scipy.optimize import minimize


r_g_sun = ((2 * G * M_sun) / (c**2)) / u.meter  # [m] Schwarzschild radius of the sun
F0_sun = (R_sun / u.meter)**2 / (2 * r_g_sun)  # [m] closest focus of SGL, ~547.303 au
v_critical = 122.3  # [GHz] Returns the frequency at which the sun has no focus


def d_critical(frequency):
    """Returns focal distance as function of frequency [GHz],
    for gravity + plasma"""

    def distance_min(impact_parameter):
        return F0_sun * impact_parameter**2 * (1 - (v_critical**2 /
            frequency ** 2) * (1 / impact_parameter)**15)**-1

    return minimize(distance_min, bounds=[(1, 2)], x0=1.1).fun[0]


def holevo_perfect(M, eta):
    """Holevo capacity of a noise free channel.
    M: # Number of photons per mode (the signal)
    eta: Receiver efficiency"""

    return ((1 + M * eta) * log2(1 + M * eta) - M * eta * log2(M * eta)) / M * eta


def holevo_thermal(M, N_th, eta):
    """Holevo capacity of a channel with noise. eta: Receiver efficiency.
    Example:
    F = 100  # number of photons
    N = 100  # Number of noise photons
    modes = 10**5  # number of modes
    M = F / modes  # Number of photons per mode (the signal)
    N_th = N / modes  # average number of noise photons per mode (the noise)"""

    def gx(x, eta, M):
        return (((1 + x)) * log2(1 + x) - x * log2(x)) / M * eta

    a = M * eta + (1 - eta) * N_th
    b = (1 - eta) * N_th
    return gx(a, eta, M) - gx(b, eta, M)


def photons_received(D_r, D_t, P, wavelength, R, Q_R=1.22):
    """Number of photons that telescope with aperture D_r [m] receives,
    D_t [m] aperture of transmitting telescope,
    wavelength lambda in [m],
    R [m] distance between D_r and D_t
    Q_R [1] diffraction limit"""

    # wavelength = wavelength * 10**-9  # nm to m
    h = 6.62607004E-34  # [J*s] Planck's h
    # c = 299792458  # [m/s] speed of light
    f = 1 / (wavelength / c) / (u.meter / u.second)
    # print(f)
    F_r = P / ((pi * h * f) * (Q_R * wavelength / D_t * R)**2) * pi * D_r**2 / 4
    return F_r


def Q(wavelength):
    """Returns the technologically achievable focusing (Q-factor)
    (Earth 2007 tech) for a given wavelength (in m)"""
    wavelength = wavelength / 10**-9
    if wavelength > 300:
        Q = 1.22
    else:  # power-law fit to values from Earth 2017 technology (see Figure)
        Q = 1 / ((5.21575598 * (wavelength / 1000)**1.20669934) / 1.22) + 0.22
    return Q


def sgl_psf(rho, wavelength, z, r_g=r_g_sun):
    """Returns the gain for a point of the PSF (rho: distance from axis)"""
    return 4 * pi**2 / (1-exp(1) ** (-4 * pi**2 * r_g / wavelength)) * r_g / \
        wavelength * j0(2*pi*(rho/wavelength)*sqrt((2*r_g)/z))**2


def integrated_flux(wavelength, z, d0, r_g=r_g_sun):
    """Eq. 143 in Turyshev (2017)"""
    first_term = 4 * pi**2 / (1 - exp(1) ** (-4 * pi**2 * r_g / wavelength)) * \
        r_g / wavelength
    zeroth_bessel = j0(pi*(d0/wavelength)*sqrt((2*r_g)/z))**2
    first_bessel =  j1(pi*(d0/wavelength)*sqrt((2*r_g)/z))**2
    return first_term * (zeroth_bessel + first_bessel)


def classical_aperture(wavelength, z, d0, r_g=r_g_sun):
    """Returns the corresponding size D of a classical aperture [m] for a
    given SGL aperture d0 [m]"""
    flux = integrated_flux(wavelength, z, d0, r_g=r_g)
    sgl_aperture = pi * (d0 / 2)**2
    effective_flux = flux * sgl_aperture
    classical_aperture = 2 * sqrt(effective_flux / pi)  # equivalent telescope diameter [m]
    return classical_aperture


def b_of_z(z, F0=F0_sun):
    """Returns impact parameter b as a function of heliocentric distance z [m]"""
    return sqrt(z) / sqrt(F0)


def z_of_b(b, F0=F0_sun):
    """Returns heliocentric distance z [m] as a function of impact parameter b"""
    return F0 * (b * (R_sun / u.meter))**2 / (R_sun / u.meter)**2


def noise_photons(wavelength, z, d0, fractional_bandwidth, print_status=False):
    """Noise photons from ring-shaped area in solar corona overlaying the
    Einstein ring"""

    def apparent_size_sun(z):
        """Returns the apparent diameter of the sun depending on
        heliocentric distance z """
        return 3600 * 180 / pi * 2 * R_sun / u.meter / z  # arcseconds

    def angular_resolution(d, wavelength, Q):
        """Returns the angular resolution of a telescope
        d: circular aperture diameter [m]
        wavelength \lambda [m]
        Q: Quality < 1.22, diffraction-limit: Q=1.22"""
        return 3600 * 180 / pi * Q * wavelength / d

    b = b_of_z(z=z)
    width_sun = apparent_size_sun(z=z)
    resolution = angular_resolution(d=d0, wavelength=wavelength, Q=1.22)
    r = b * width_sun / 2
    A_ring = pi * ((r + resolution)**2 - r**2)
    A_sun = pi * width_sun**2 / 4
    ring_brightness = b**-6 * 10**-6
    noise_in_ring = ring_brightness * (A_ring / A_sun)  # [solar luminosities]
    distance_scale = 4 * pi * z**2
    noise_in_ring_scaled_to_distance = noise_in_ring / distance_scale
    noise_watt = noise_in_ring_scaled_to_distance * L_sun / u.watt
    joule_wavelength = 1.98641E-19 / wavelength / 10**6
    fractional_noise = fractional_bandwidth * noise_watt
    noise_photons_per_nm_per_sec = fractional_noise / joule_wavelength

    if print_status:
        print('b=', b)
        print('width of sun (arcsec)', width_sun)
        print('resolution', resolution)
        print('r', r)
        print('A_ring', A_ring)
        print('A_sun', A_sun)
        print('A_ring / A_sun', A_ring / A_sun)
        print('ring_brightness', ring_brightness)
        print('noise_in_ring_scaled_to_distance', noise_in_ring_scaled_to_distance)
        print('noise_watt', noise_watt)
        print('joule_wavelength', joule_wavelength)
        print('fractional_bandwidth', fractional_bandwidth)
        print('fractional_noise', fractional_noise)
        print('noise_photons_per_nm_per_sec', noise_photons_per_nm_per_sec)

    return noise_photons_per_nm_per_sec


def QuadraticLimbDarkening(Impact, limb1, limb2):
    """Quadratic limb darkening. Kopal 1950, Harvard Col. Obs. Circ., 454, 1"""
    return 1 - limb1 * (1 - Impact) - limb2 * (1 - Impact) ** 2


def xy_rotate(x, y, xcen, ycen, phi):
    """
    Copyright 2009 by Adam S. Bolton
    Creative Commons Attribution-Noncommercial-ShareAlike 3.0 license applies
    NAME: xy_rotate

    PURPOSE: Transform input (x, y) coordiantes into the frame of a new
             (x, y) coordinate system that has its origin at the point
             (xcen, ycen) in the old system, and whose x-axis is rotated
             c.c.w. by phi degrees with respect to the original x axis.

    USAGE: (xnew,ynew) = xy_rotate(x, y, xcen, ycen, phi)

    ARGUMENTS:
      x, y: numpy ndarrays with (hopefully) matching sizes
            giving coordinates in the old system
      xcen: old-system x coordinate of the new origin
      ycen: old-system y coordinate of the new origin
      phi: angle c.c.w. in degrees from old x to new x axis

    RETURNS: 2-item tuple containing new x and y coordinate arrays

    WRITTEN: Adam S. Bolton, U. of Utah, 2009
    """
    phirad = numpy.deg2rad(phi)
    xnew = (x - xcen) * numpy.cos(phirad) + (y - ycen) * numpy.sin(phirad)
    ynew = (y - ycen) * numpy.cos(phirad) - (x - xcen) * numpy.sin(phirad)
    return (xnew,ynew)


def gauss_2d(x, y, par):
    """
    Copyright 2009 by Adam S. Bolton
    Creative Commons Attribution-Noncommercial-ShareAlike 3.0 license applies
    NAME: gauss_2d

    PURPOSE: Implement 2D Gaussian function

    USAGE: z = gauss_2d(x, y, par)

    ARGUMENTS:
      x, y: vecors or images of coordinates;
            should be matching numpy ndarrays
      par: vector of parameters, defined as follows:
        par[0]: amplitude
        par[1]: intermediate-axis sigma
        par[2]: x-center
        par[3]: y-center
        par[4]: axis ratio
        par[5]: c.c.w. major-axis rotation w.r.t. x-axis

    RETURNS: 2D Gaussian evaluated at x-y coords

    NOTE: amplitude = 1 is not normalized, but rather has max = 1

    WRITTEN: Adam S. Bolton, U. of Utah, 2009
    """
    (xnew,ynew) = xy_rotate(x, y, par[2], par[3], par[5])
    r_ell_sq = ((xnew**2)*par[4] + (ynew**2)/par[4]) / numpy.abs(par[1])**2
    return par[0] * numpy.exp(-0.5*r_ell_sq)


def sie_grad(x, y, par):
    """
    Copyright 2009 by Adam S. Bolton
    Creative Commons Attribution-Noncommercial-ShareAlike 3.0 license applies
    NAME: sie_grad

    PURPOSE: compute the deflection of an SIE potential

    USAGE: (xg, yg) = sie_grad(x, y, par)

    ARGUMENTS:
      x, y: vectors or images of coordinates;
            should be matching numpy ndarrays
      par: vector of parameters with 1 to 5 elements, defined as follows:
        par[0]: lens strength, or 'Einstein radius'
        par[1]: (optional) x-center (default = 0.0)
        par[2]: (optional) y-center (default = 0.0)
        par[3]: (optional) axis ratio (default=1.0)
        par[4]: (optional) major axis Position Angle
                in degrees c.c.w. of x axis. (default = 0.0)

    RETURNS: tuple (xg, yg) of gradients at the positions (x, y)

    NOTES: This routine implements an 'intermediate-axis' conventionumpy.
      Analytic forms for the SIE potential can be found in:
        Kassiola & Kovner 1993, ApJ, 417, 450
        Kormann et al. 1994, A&A, 284, 285
        Keeton & Kochanek 1998, ApJ, 495, 157
      The parameter-order convention in this routine differs from that
      of a previous IDL routine of the same name by ASB.

    WRITTEN: Adam S. Bolton, U of Utah, 2009
    """
    # Set parameters:
    b = numpy.abs(par[0]) # can't be negative!!!
    xzero = 0. if (len(par) < 2) else par[1]
    yzero = 0. if (len(par) < 3) else par[2]
    q = 1. if (len(par) < 4) else numpy.abs(par[3])
    phiq = 0. if (len(par) < 5) else par[4]
    eps = 0.001 # for sqrt(1/q - q) < eps, a limit expression is used.
    # Handle q > 1 gracefully:
    if (q > 1.):
        q = 1.0 / q
        phiq = phiq + 90.0
    # Go into shifted coordinats of the potential:
    phirad = numpy.deg2rad(phiq)
    xsie = (x-xzero) * numpy.cos(phirad) + (y-yzero) * numpy.sin(phirad)
    ysie = (y-yzero) * numpy.cos(phirad) - (x-xzero) * numpy.sin(phirad)
    # Compute potential gradient in the transformed system:
    r_ell = numpy.sqrt(q * xsie**2 + ysie**2 / q)
    qfact = numpy.sqrt(1./q - q)
    # (r_ell == 0) terms prevent divide-by-zero problems
    if (qfact >= eps):
        xtg = (b/qfact) * numpy.arctan(qfact * xsie / (r_ell + (r_ell == 0)))
        ytg = (b/qfact) * numpy.arctanh(qfact * ysie / (r_ell + (r_ell == 0)))
    else:
        xtg = b * xsie / (r_ell + (r_ell == 0))
        ytg = b * ysie / (r_ell + (r_ell == 0))
    # Transform back to un-rotated system:
    xg = xtg * numpy.cos(phirad) - ytg * numpy.sin(phirad)
    yg = ytg * numpy.cos(phirad) + xtg * numpy.sin(phirad)
    # Return value:
    return (xg, yg)
