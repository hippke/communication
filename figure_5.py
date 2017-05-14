import numpy
import matplotlib.pyplot as plt
import matplotlib.patches as patches


path = ''
data = numpy.genfromtxt(
    path + 'data_atmo_glow.txt',
    dtype=[
        ('wavelength', 'f8'),
        ('glow', 'f8')])

plt.rc('font',  family='serif', serif='Computer Modern Roman')
plt.rc('text', usetex=True)
size = 4.5
aspect_ratio = 1.5
plt.figure(figsize=(size, size / aspect_ratio))
ax = plt.gca()
ax.get_yaxis().set_tick_params(direction='out')
ax.get_xaxis().set_tick_params(direction='out')
ax.get_yaxis().set_tick_params(which='both', direction='out')
ax.get_xaxis().set_tick_params(which='both', direction='out')
ax.set_xlabel(r'Wavelength $\lambda$ ($\mu$m)')
ax.set_ylabel(r'Flux ($\gamma$\,nm$^{-1}$\,s$^{-1}$\,arcsec$^{-2}$\,m$^{-2}$)')
ax.set_xlim(0.3, 30)
ax.set_ylim(0.001, 20000000)
ax.set_xscale('log')
ax.set_yscale('log')
plt.plot(data['wavelength'], data['glow'], color='black', linewidth=0.25)

cm = plt.cm.get_cmap('jet')
for i in range(1000):
    ax.add_patch(patches.Rectangle((0.35 + i/1000/3, 1e-3), 0.01, 0.005,
        facecolor=cm(i/1000), edgecolor='none'))
ax.text(0.72, 1.5e-3, r'IR')

plt.savefig(path + 'figure_atmo_glow.pdf', bbox_inches='tight')
#plt.show()
