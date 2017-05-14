import numpy
import matplotlib.pyplot as plt
import matplotlib.patches as patches


path = ''
data = numpy.genfromtxt(
    path + 'data_atmo.csv',
    dtype=[
        ('wavelength', 'f8'),
        ('transmission', 'f8')])

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
ax.set_ylabel(r'Atmospheric transmission')
ax.set_xlim(0.1, 20000)
ax.set_ylim(0, 1.1)
ax.set_xscale('log')
plt.plot(data['wavelength'], data['transmission'], color='black', linewidth=0.5)

cm = plt.cm.get_cmap('jet')
for i in range(1000):
    ax.add_patch(patches.Rectangle((0.35 + i/1000/3, 1), 0.01, 0.1,
        facecolor=cm(i/1000), edgecolor='none'))
ax.text(0.12, 1.02, r'UV')
ax.text(0.80, 1.02, r'IR')

plt.savefig(path + 'figure_atmo.pdf', bbox_inches='tight')

data = numpy.genfromtxt(
    path + 'data_transparency.txt',
    dtype=[
        ('wavelength', 'f8'),
        ('transmission', 'f8')])

plt.clf()
plt.figure(figsize=(size, size / aspect_ratio))
ax = plt.gca()
ax.get_yaxis().set_tick_params(direction='out')
ax.get_xaxis().set_tick_params(direction='out')
ax.get_yaxis().set_tick_params(which='both', direction='out')
ax.get_xaxis().set_tick_params(which='both', direction='out')
ax.set_xlabel(r'Wavelength $\lambda$ (nm)')
ax.set_ylabel(r'Atmospheric transmission')
ax.set_xlim(930, 940)
ax.set_ylim(0, 1.1)
plt.plot(data['wavelength'] * 1000, data['transmission'], color='black')

plt.savefig(path + 'figure_atmo_zoom.pdf', bbox_inches='tight')
#plt.show()
