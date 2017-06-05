import numpy
import matplotlib.patches as patches
from scipy.optimize import minimize
from astropy.constants import c, G, M_sun, R_sun, au
from matplotlib import pyplot as plt


def d_critical(frequency):
    """Returns focal distance as function of frequency [GHz],
       for gravity + plasma"""

    def distance_min(impact_parameter):
        return d_focal * impact_parameter**2 * (1 - (v_critical**2 /
            frequency**2) * (1 / impact_parameter)**15)**-1

    return minimize(distance_min, bounds=[(1, 2)], x0=1.1).fun[0]



v_critical = 122.3  # GHz
schwarzschild_radius_sun = (2 * G * M_sun) / (c**2)
d_focal = (R_sun**2) / (2 *  schwarzschild_radius_sun) / au
gridsize = 1000
freqs = numpy.logspace(numpy.log10(v_critical), 4, gridsize, endpoint=True)

#print(freqs)
distances = numpy.zeros(gridsize)
print(d_focal)


#for i in range(len())


#print((1.1*R_sun**2) / (2 *  schwarzschild_radius_sun) / au)
#print(schwarzschild_radius_sun)

counter = 0
for freq in freqs:
    distances[counter] = d_critical(freq)
    counter += 1


plt.rc('font',  family='serif', serif='Computer Modern Roman')
plt.rc('text', usetex=True)
size = 4
aspect_ratio = 1.5
plt.figure(figsize=(size, size / aspect_ratio))
ax = plt.gca()
fig, ax = plt.subplots(figsize=(size, size / aspect_ratio))
ax2 = ax.twinx()
ax.get_yaxis().get_major_formatter().set_useOffset(False)
ax.get_yaxis().get_major_formatter().set_scientific(False)
ax.get_yaxis().set_tick_params(direction='out')
ax.get_xaxis().set_tick_params(direction='out')
ax.get_yaxis().set_tick_params(which='both', direction='out')
ax.get_xaxis().set_tick_params(which='both', direction='out')
ax2.get_yaxis().get_major_formatter().set_useOffset(False)
ax2.get_yaxis().get_major_formatter().set_scientific(False)
ax2.get_yaxis().set_tick_params(direction='out')
ax2.get_xaxis().set_tick_params(direction='out')
ax2.get_yaxis().set_tick_params(which='both', direction='out')
ax2.get_xaxis().set_tick_params(which='both', direction='out')

ax2.set_ylim([0.9613, 1.24622])
ax2.set_ylabel(r'Distance of Einstein ring ($b/R_{\odot}$)')

ax.set_xscale('log')
ax.add_patch(patches.Rectangle((100, 500), v_critical - 100, 500,  fill=False, hatch='//'))
ax.plot(freqs, distances, color='black')
ax.plot((100, 10**4), (d_focal, d_focal), color='black', linestyle = 'dashed')
ax.scatter(410, 600, s=100, marker='*')
ax.text(170, 552, '547.3 au')
ax.text(125, 505, '122.3 GHz')
ax.set_xlabel('Frequency (GHz)')
ax.set_ylabel('Distance from Sun (au)')
ax.set_ylim([500, 850])
ax.set_xlim([100, 10000])

plt.savefig('figure_distance_a.pdf', bbox_inches='tight')

plt.show()
