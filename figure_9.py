from PyCom import holevo_perfect, holevo_thermal
import numpy
from math import log2
from matplotlib import pyplot as plt


path = ''
plt.rc('font',  family='serif', serif='Computer Modern Roman')
plt.rc('text', usetex=True)
size = 4.5
aspect_ratio = 1.5
plt.figure(figsize=(size, size / aspect_ratio))
ax = plt.gca()

eta = 0.5 # channel power transmissivity
N = 10**-5
grid_size = 10000
grid = numpy.logspace(-6, 3, num=grid_size)
grid_values_thermal = numpy.zeros([grid_size])

for value in range(len(grid)):
    grid_values_thermal[value] = holevo_thermal(N=N, N_th=grid[value], eta=eta)
ax.plot(grid_values_thermal, grid, color='black')

# Intersection
ax.plot([1, 1], [10**-6, 0.133], color='black', linestyle='dashed')
ax.plot([10**-4, 1], [0.133,0.133], color='black', linestyle='dashed')
ax.get_yaxis().get_major_formatter().set_useOffset(False)
ax.get_yaxis().get_major_formatter().set_scientific(False)
ax.get_yaxis().set_tick_params(direction='out')
ax.get_xaxis().set_tick_params(direction='out')
ax.get_yaxis().set_tick_params(which='both', direction='out')
ax.get_xaxis().set_tick_params(which='both', direction='out')
ax.set_xscale('log')
ax.set_yscale('log')
ax.text(0.8, 2.25 * 10**-2, r'$C_{th} = 1$\,bits/photon possible for', color='black', horizontalalignment='right')
ax.text(0.8, 2.25 * 10**-3, r'$N_{M} = 0.13$\,noise photons/mode', color='black', horizontalalignment='right')
ax.text(0.8, 2.25 * 10**-4, r'with $M=10^{-5}$, $\eta=0.5$', color='black', horizontalalignment='right')
ax.set_xlabel('Capacity (bits/photon)')
ax.set_ylabel('Thermal noise $N_M$ (photons/mode)')
plt.savefig(path + 'figure_holevo_noise.pdf', bbox_inches='tight')
#plt.show()
