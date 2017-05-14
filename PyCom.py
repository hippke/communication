from math import log, log2, pi


def holevo_perfect(N, eta):
    return ((1 + N * eta) * log2(1 + N * eta) - N * eta * log2(N * eta)) / N * eta


def gx(x, eta, N):
    return (((1 + x)) * log2(1 + x) - x * log2(x)) / N * eta


def holevo_thermal(N, N_th, eta):
    a = N * eta + (1 - eta) * N_th
    b = (1 - eta) * N_th
    return gx(a, eta, N) - gx(b, eta, N)


def photons_received(D_r, D_t, P, wavelength, R, Q_R=1.22):
    # \lambda refered to as "wavelength" because it is a protected word in Python
    wavelength = wavelength * 10**-9  # nm to m
    h = 6.62607004E-34  # Js
    c = 299792458  # m/s
    f = 1 / (wavelength / c)
    F_r = P / ((pi * h * f) * (Q_R * wavelength / D_t * R)**2) * pi * D_r**2 / 4
    return F_r


def Q(wavelength):
    """Returns the technologically achievable focusing (Q-factor)
    (Earth 2007 tech) for a given wavelength (in nm)"""
    if wavelength > 300:
        Q = 1.22
    else:  # power-law fit to values from Earth 2017 technology (see Figure)
        Q = 1 / ((5.21575598 * (wavelength / 1000)**1.20669934) / 1.22) + 0.22
    return Q

