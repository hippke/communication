import numpy
import matplotlib.patches as patches
from scipy.optimize import minimize
from astropy.constants import c, G, M_sun, R_sun, au
from matplotlib import pyplot as plt


def QuadraticLimbDarkening(Impact, limb1, limb2):
    """Quadratic limb darkening. Kopal 1950, Harvard Col. Obs. Circ., 454, 1"""
    return 1 - limb1 * (1 - Impact) - limb2 * (1 - Impact) ** 2


path = '' # 'I:/D/Research/ET sail/MyPaper/part2/'
plt.rc('font',  family='serif', serif='Computer Modern Roman')
plt.rc('text', usetex=True)
size = 11
aspect_ratio = 1.5
plt.figure(figsize=(size, size / aspect_ratio))

ax = plt.subplot(121, aspect='equal')


# Left plot
#plt.subplot(121, aspect='equal')
ax = plt.gca()
ax.set_ylim([-1.6, 1.6])
ax.set_xlim([-1.6, 1.6])
ax.get_yaxis().get_major_formatter().set_useOffset(False)
ax.get_yaxis().get_major_formatter().set_scientific(False)
ax.get_yaxis().set_tick_params(direction='out')
ax.get_xaxis().set_tick_params(direction='out')
ax.get_yaxis().set_tick_params(which='both', direction='out')
ax.get_xaxis().set_tick_params(which='both', direction='out')
ax.tick_params(direction='out')
ax.tick_params(axis='both', which='major')
background = plt.Circle((0, 0), 3, color='black', fill=True)
plt.gcf().gca().add_artist(background)


# Corona
Quality = 100  # Radius of star in pixels --> size of numerical sampling grid
for i in range(Quality):
    c = float(i/100.)
    corona = plt.Circle((0, 0), 2 - i / float(Quality),
        color=(c, c, c))
    plt.gcf().gca().add_artist(corona)

# Sun
limb1 = 0.5971
limb2 = 0.1172
for i in range(Quality):
    Impact = i / float(Quality)
    LimbDarkening = QuadraticLimbDarkening(Impact, limb1, limb2)
    Sun = plt.Circle((0, 0), 1 - i / float(Quality),
        color=(LimbDarkening, LimbDarkening, 0))
    plt.gcf().gca().add_artist(Sun)

einstein_ring = plt.Circle((0, 0), 1.35, color='white', fill=False)
plt.gcf().gca().add_artist(einstein_ring)

#ax.text(-1.55, 1.4, r'$z=1000$\,au', color='white')
#ax.text(-1.55, 1.25, r'$b=1.35$', color='white')

#ax.text(0.6, 1.30, 'Einstein ring', color='white')
#ax.set_xlabel('Distance (stellar radii)', color='white')
#ax.set_ylabel('Distance (stellar radii)', color='white')
ax.set_xlabel('Distance (stellar radii)')
ax.set_ylabel('Distance (stellar radii)')
# ax.arrow(0, 0, 0, 1.30, head_width=0.05, head_length=0.05, color='green')


# Right plot
#plt.subplot(121, aspect='equal')
ax2 = plt.subplot(122, aspect='equal')
#ax2.set_ylim([0, 1.5])
#ax2.set_xlim([0, 1.5])
ax2.set_xlabel('Distance (stellar radii)')
#ax2.set_ylabel('Distance (stellar radii)')
ax2.set_ylim([-1.6, 1.6])
ax2.set_xlim([-1.6, 1.6])
ax2.get_yaxis().get_major_formatter().set_useOffset(False)
ax2.get_yaxis().get_major_formatter().set_scientific(False)
ax2.get_yaxis().set_tick_params(direction='out')
ax2.get_xaxis().set_tick_params(direction='out')
ax2.get_yaxis().set_tick_params(which='both', direction='out')
ax2.get_xaxis().set_tick_params(which='both', direction='out')
ax2.tick_params(direction='out')
ax2.tick_params(axis='both', which='major')
background = plt.Circle((0, 0), 3, color='black', fill=True)
plt.gcf().gca().add_artist(background)

"""
# outer_coronograph
for i in range(Quality):
    c = float(i/100.)
    outer_coronograph = plt.Circle((0, 0), 1.3 + i / float(Quality),
        color=(0, 0, 0))
plt.gcf().gca().add_artist(outer_coronograph)
"""

# Corona
Quality = 100  # Radius of star in pixels --> size of numerical sampling grid
for i in range(Quality):
    if i < 50:
        c = 0.  # black, outer coronograph
    else:
        c = float(i/100.)
    corona = plt.Circle((0, 0), 2 - i / float(Quality),
        color=(c, c, c))
    plt.gcf().gca().add_artist(corona)



coronograph = plt.Circle((0, 0), 1.2, color='black', fill=True)
plt.gcf().gca().add_artist(coronograph)

einstein_ring = plt.Circle((0, 0), 1.35, color='white', fill=False)
plt.gcf().gca().add_artist(einstein_ring)

#aperture_ring = plt.Circle((-1.25, -1.18), 0.0588 / 2, color='white', fill=True)
#plt.gcf().gca().add_artist(aperture_ring)
#plt.subplot(122, aspect='equal')
#plt.plot([1, 2, 3], [1, 2, 3])
#ax2.text(-1.15, -1.4, r'Airy disk for $d=1$m, $\lambda=1\mu$m,' '\n' r'$\theta=0.175$ arcsec', color='white')


plt.savefig(path + 'sun_double.pdf', bbox_inches='tight')
