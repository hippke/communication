import numpy
import matplotlib.pyplot as plt
import matplotlib.patches as patches


path = ''
data_cmb = numpy.genfromtxt('data_cmb.txt', dtype=[('wavelength', 'f8'),('flux', 'f8')])
data_ir = numpy.genfromtxt('data_ir.txt', dtype=[('wavelength', 'f8'),('flux', 'f8')])
data_xray = numpy.genfromtxt('data_xray.txt', dtype=[('wavelength', 'f8'),('flux', 'f8')])
data_gammaray = numpy.genfromtxt('data_gammaray.txt', dtype=[('wavelength', 'f8'),('flux', 'f8')])

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
ax.set_xlim(5e-7, 1e5)
ax.set_ylim(4e-10, 2e-4)
ax.set_xscale('log')
ax.set_yscale('log')

unit_conversion = 1.99E-07
# from nW\,m$^{-2}$\,sr${^-2}$ to $\gamma$\,nm$^{-1}$\,s$^{-1}$\,arcsec$^2$\,m$^{-2}$

plt.plot(data_cmb['wavelength'][55:], data_cmb['flux'][55:] * unit_conversion, color='black')
plt.plot(data_ir['wavelength'][:178], data_ir['flux'][:178] * unit_conversion, color='blue')
plt.scatter(data_xray['wavelength'], data_xray['flux'] * unit_conversion, color='orange')
plt.scatter(data_gammaray['wavelength'], data_gammaray['flux'] * unit_conversion, color='red')

ax.text(1e-6, 1 * unit_conversion, r'$\gamma$-ray', color='red')
ax.text(1e-4, 1 * unit_conversion, r'X-ray', color='orange')
ax.text(1e+3, 1 * unit_conversion, r'CMB', color='black')

cm = plt.cm.get_cmap('jet')
for i in range(1000):
    ax.add_patch(patches.Rectangle((0.35 + i/1000/3, 3e-3 * unit_conversion), 0.01, 0.006 * unit_conversion,
        facecolor=cm(i/1000), edgecolor='none'))
ax.text(0.045, 4e-3 * unit_conversion, r'UV')
ax.text(0.90, 4e-3 * unit_conversion, r'IR')

plt.savefig(path + 'figure_background_sources.pdf', bbox_inches='tight')
#plt.show()
