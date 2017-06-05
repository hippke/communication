import numpy
import matplotlib.pyplot as plt
import matplotlib.patches as patches


path = ''
data_all = numpy.genfromtxt('data_AllMODEtr.txt',
    usecols=(1,3),
    dtype=[('wavelength', 'f8'),('flux', 'f8')])

data_uv = numpy.genfromtxt('data_uv.txt',
    usecols=(2,6),
    dtype=[('wavelength', 'f8'),('flux', 'f8')])

data_bass_uv = numpy.genfromtxt('data_bass_uv.txt',
    dtype=[('wavelength', 'f8'),('flux', 'f8')])

data_bass2000 = numpy.genfromtxt('data_bass2000.txt',
    dtype=[('wavelength', 'f8'),('flux', 'f8')])


solar_constant = 1366.91
max_peak = numpy.nanmax(data_all['flux'])
solar_luminosity = 3.828 * 10**26  # Watt
scale_to_fraction = max_peak / solar_constant * solar_luminosity



# print(max_peak, scale_to_fraction)

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
ax.set_xlabel(r'Wavelength $\lambda$ (nm)')
ax.set_ylabel(r'Solar luminosity (W\,nm$^{-1}$)')
ax.set_xlim(100, 1 * 10**5)
ax.set_ylim(1e16, 2e24)
ax.set_xscale('log')
ax.set_yscale('log')

ax.text(300, 1e20, r'$L_{\odot}=3.828\times10^{26}$\,W')
ax.text(14000, 5*1e23, r'$\lambda>12$\,$\mu$m')
ax.text(14000, 7*1e22, 'starshade')
ax.text(14000, 9*1e21, 'required for')
ax.text(14000, 11*1e20, r'$d<1$\,m')


uv_scale = 4.0041e-05 / 100  # they nornalized each UV segment to 1.0

plt.plot((12000, 12000), (10**16, 10**26), color='black', linestyle='dashed')
#plt.plot((120000, 120000), (10**17, 10**19), color='blue')

#plt.plot(data_all['wavelength'][:49169], data_all['flux'][:49169] * scale_to_fraction, color='black', linewidth=0.5)
plt.plot(data_all['wavelength'], data_all['flux'] * scale_to_fraction, color='black', linewidth=0.5)
plt.plot(data_uv['wavelength'], data_uv['flux'] * scale_to_fraction, color='black', linewidth=0.5)
plt.plot(data_bass_uv['wavelength']/10, data_bass_uv['flux'] * uv_scale * scale_to_fraction, color='black', linewidth=0.5)


#plt.plot(data_all['wavelength'][49169:], data_all['flux'][49169:] *10**6 * scale_to_fraction, color='black', linewidth=0.5)


cm = plt.cm.get_cmap('jet')
for i in range(1000):
    ax.add_patch(patches.Rectangle((0.35*1000 + i/3, 10**16), 1, 4 * 10**16,
        facecolor=cm(i/1000), edgecolor='none'))
ax.text(200, 1.5 * 10**16, r'UV')
ax.text(750, 1.5 * 10**16, r'IR')


plt.savefig(path + 'figure_spectrum_all.pdf', bbox_inches='tight')
#plt.show()
"""
plt.clf()
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
ax.set_xlabel(r'Wavelength $\lambda$ (nm)')
ax.set_ylabel(r'Irradiance (W\,m$^{-2}$\,nm$^{-1}$)')
ax.set_xlim(654.4, 657.9)
ax.set_ylim(0, 1.1*10**24)
ax.set_xlabel(r'Wavelength $\lambda$ (nm)')
ax.set_ylabel(r'Irradiance (W\,nm$^{-1}$)')
scale = 6434


ax.text(654.5, 0.05*10**24, r'Ti I')
ax.text(656.15, 0.05*10**24, r'H I')
ax.text(656.8, 0.05*10**24, r'Sc I')
ax.text(657.4, 0.05*10**24, r'Fe I')

ax.get_yaxis().get_major_formatter().set_useOffset(False)
#ax.get_yaxis().get_major_formatter().set_scientific(False)
plt.plot(data_bass2000['wavelength']/10, data_bass2000['flux'] / scale * scale_to_fraction, color='black', linewidth=1)
plt.savefig(path + 'figure_spectrum_6562.pdf', bbox_inches='tight')
"""