from PyCom import holevo_perfect, holevo_thermal
import numpy
from matplotlib import pyplot as plt


plt.rc('font',  family='serif', serif='Computer Modern Roman')
plt.rc('text', usetex=True)
size = 4.5
aspect_ratio = 1.5
plt.figure(figsize=(size, size / aspect_ratio))
ax = plt.gca()

eta = 0.5 # channel power transmissivity
N_th = 0.001  # noise photons per mode

grid_size = 10000
grid = numpy.logspace(-6, 1, num=grid_size)
grid_values_perfect = numpy.zeros([grid_size])
grid_values_thermal = numpy.zeros([grid_size])


for value in range(len(grid)):
    grid_values_perfect[value] = holevo_perfect(N=grid[value], eta=eta)
ax.plot(grid_values_perfect, grid, color='black')

for value in range(len(grid)):
    grid_values_thermal[value] = holevo_thermal(N=grid[value], N_th=0.01, eta=eta)
ax.plot(grid_values_thermal, grid, color='red')

for value in range(len(grid)):
    grid_values_thermal[value] = holevo_thermal(N=grid[value], N_th=0.001, eta=eta)
ax.plot(grid_values_thermal, grid, color='red')

for value in range(len(grid)):
    grid_values_thermal[value] = holevo_thermal(N=grid[value], N_th=0.0001, eta=eta)
ax.plot(grid_values_thermal, grid, color='red')


ax.get_yaxis().get_major_formatter().set_useOffset(False)
ax.get_yaxis().get_major_formatter().set_scientific(False)
ax.get_yaxis().set_tick_params(direction='out')
ax.get_xaxis().set_tick_params(direction='out')
ax.get_yaxis().set_tick_params(which='both', direction='out')
ax.get_xaxis().set_tick_params(which='both', direction='out')
ax.set_yscale('log')

ax.text(1.12, 2.5 * 10**-6, r'$N_{M}=$', color='black')
ax.text(2.00, 3 * 10**-6, r'$0.1$', color='black')
ax.text(2.82, 3 * 10**-6, r'$0.01$', color='black')
ax.text(3.65, 3 * 10**-6, r'$0.001$', color='black')
ax.text(4.00, 1 * 10**-4, r'$C_{ult}$ for $\eta=0.5$', color='black')

ax.set_ylim([1e-6, 1])
ax.set_xlim([0, 6])
ax.set_xlabel('Capacity $C_{th}$ (bits/photon)')
ax.set_ylabel('$M$ (photons/mode)')
plt.savefig('figure_holevo.pdf', bbox_inches='tight')

#plt.show()
