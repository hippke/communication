import numpy
import matplotlib.patches as patches
from scipy.optimize import minimize
from astropy.constants import c, G, M_sun, R_sun, au
from matplotlib import pyplot as plt
from PyCom import QuadraticLimbDarkening


plt.rc('font',  family='serif', serif='Computer Modern Roman')
plt.rc('text', usetex=True)
size = 5
aspect_ratio = 1.5
plt.figure(figsize=(size, size / aspect_ratio))
ax = plt.gca()
ax.get_yaxis().get_major_formatter().set_useOffset(False)
ax.get_yaxis().get_major_formatter().set_scientific(False)
ax.get_yaxis().set_tick_params(direction='out')
ax.get_xaxis().set_tick_params(direction='out')
ax.get_yaxis().set_tick_params(which='both', direction='out')
ax.get_xaxis().set_tick_params(which='both', direction='out')

limb1 = 0.5971
limb2 = 0.1172
Quality = 100  # Radius of star in pixels --> size of numerical sampling grid

for i in range(Quality):
    Impact = i / float(Quality)
    LimbDarkening = QuadraticLimbDarkening(Impact, limb1, limb2)
    Sun = plt.Circle((0, 0), 1 - i / float(Quality),
        color=(LimbDarkening, LimbDarkening, 0))
    plt.gcf().gca().add_artist(Sun)

einstein_ring = plt.Circle((-3, 0), 1.2, color='#8bb1ed', fill=False, linewidth=2.5)
einstein_ring_outer = plt.Circle((-3, 0), 1.15, color='black', fill=False, linewidth=1)
einstein_ring_inner = plt.Circle((-3, 0), 1.25, color='black', fill=False, linewidth=1)
plt.gcf().gca().add_artist(einstein_ring)
plt.gcf().gca().add_artist(einstein_ring_outer)
plt.gcf().gca().add_artist(einstein_ring_inner)

plt.plot((-3, -2), (1.15, 1.15), color='black', linewidth=0.5, linestyle='solid')
plt.plot((-3, -2), (1.25, 1.25), color='black', linewidth=0.5, linestyle='solid')

plt.plot((-0.5, 0.25), (1.15, 1.15), color='black', linewidth=0.5, linestyle='solid')
plt.plot((-0.5, 0.25), (1.25, 1.25), color='black', linewidth=0.5, linestyle='solid')

plt.plot((-3, 0.25), (-1.15, -1.15), color='black', linewidth=0.5, linestyle='solid')
plt.plot((-3, 0.25), (-1.25, -1.25), color='black', linewidth=0.5, linestyle='solid')

plt.plot((0.25, 3.75), (-1.15, 0.1), color='black', linewidth=0.5, linestyle='solid')
plt.plot((0.25, 3.75), (-1.25, 0), color='black', linewidth=0.5, linestyle='solid')

plt.plot((0.25, 3.75), (1.15, -0.1), color='black', linewidth=0.5, linestyle='solid')
plt.plot((0.25, 3.75), (1.25, 0), color='black', linewidth=0.5, linestyle='solid')

# d
plt.plot((4, 4), (-0.12, 0.12), color='black', linewidth=0.5, linestyle='solid')
plt.plot((3.95, 4.05), (-0.12, -0.12), color='black', linewidth=0.5, linestyle='solid')
plt.plot((3.95, 4.05), (+0.12, +0.12), color='black', linewidth=0.5, linestyle='solid')
ax.text(4.15, -0.1, r'$d$')

# z
plt.plot((0, 3.85), (-1.65, -1.65), color='black', linewidth=0.5, linestyle='solid')
plt.plot((0, 0), (-1.6, -1.7), color='black', linewidth=0.5, linestyle='solid')
plt.plot((3.85, 3.85), (-1.6, -1.7), color='black', linewidth=0.5, linestyle='solid')
ax.text(1.75, -1.55, r'$z$')

# b
ax.arrow(-3, 0, 0, 1, fc="k", ec="k", head_width=0.05, head_length=0.05)
ax.text(-2.9, 0.3, r'$b$')

# R_*
ax.arrow(0, 0, 0, 0.85, fc="k", ec="k", head_width=0.05, head_length=0.05)
ax.text(0.1, 0.3, r'$R_{\odot}$')

plt.xticks([])
plt.yticks([])
plt.axes().set_aspect('equal')
ax.set_xlim([-5, 5])
ax.set_ylim([-2, 2])
ax.text(-4, 1.5, 'Einstein ring')
ax.text(-1.85, 0.94, r'$w=\frac{d}{2}$')
ax.tick_params(direction='out')
ax.tick_params(axis='both', which='major')
plt.savefig('sgl_scheme.pdf', bbox_inches='tight')
