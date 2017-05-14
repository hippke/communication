import numpy
import matplotlib.pyplot as plt
import matplotlib.patches as patches


path = ''
gal_center = numpy.genfromtxt(
    path + 'data_apj394440t8_mrt.txt',
    skip_header=21,
    delimiter=[8, 9, 17, 19],
    autostrip=True,
    dtype=[
        ('wavelength', 'f8'),
        ('extinction', 'f8'),
        ('uncertainty', 'f8')])

draine = numpy.genfromtxt(
    path + 'data_kext_albedo_WD_MW_3.1_60_D03.all',
    skip_header=80,
    usecols=(0,1,2,3,4,5),
    autostrip=True,
    dtype=[
        ('wavelength', 'f8'),
        ('albedo', 'f8'),
        ('cos', 'f8'),
        ('C_ext', 'f8'),
        ('K_abs', 'f8'),
        ('cos2', 'f8')])

lyman_extinct = numpy.genfromtxt(
    path + 'data_lyman_extinct.csv',
    dtype=[('factor', 'f8')])


def get_extinction(distance, wavelength):
    scale_factor = 0.2665
    A = draine['C_ext'] / scale_factor * 10**21  # normalized to 1000 pc
    curve = 1 / 2.512**(A/(1000/distance))
    curve = curve / lyman_extinct['factor']
    idx = numpy.argmin(numpy.abs(draine['wavelength'] - wavelength / 1000))
    return curve[idx]


# print(get_extinction(distance=100, wavelength=9400))

scale_factor = 0.2665
A = draine['C_ext'] / scale_factor * 10**21  # normalized to 1000 pc

distance = 10  # pc
fraction_1kpc = 1 / 2.512**A
fraction_100pc = 1 / 2.512**(A/10)
distance = 10  # pc
fraction_10pc = 1 / 2.512**(A/(1000/distance))
distance = 1  # pc
fraction_1pc = 1 / 2.512**(A/(1000/distance))
distance = 8000  # pc
fraction_1pc = fraction_1pc / lyman_extinct['factor']
fraction_10pc = fraction_10pc  / lyman_extinct['factor']
fraction_100pc = fraction_100pc / lyman_extinct['factor']
fraction_1kpc = fraction_1kpc / lyman_extinct['factor']


plt.rc('font',  family='serif', serif='Computer Modern Roman')
plt.rc('text', usetex=True)
size = 4.5
aspect_ratio = 1.5
plt.figure(figsize=(size, size / aspect_ratio))
ax = plt.gca()

ax.plot(gal_center['wavelength'], 1 / 2.512 ** gal_center['extinction'], color='black')

ax.plot(draine['wavelength'], fraction_1pc, color='green')
ax.plot(draine['wavelength'], fraction_10pc, color='orange')
ax.plot(draine['wavelength'], fraction_100pc, color='red')
ax.plot(draine['wavelength'], fraction_1kpc, color='blue')


ax.text(160, 0.04, '8kpc\ngalactic\ncenter', color='black', horizontalalignment='right')
#ax.text(32, 0.23, '8 kpc', color='black')

ax.text(0.0035, 0.40, '1 kpc', color='blue')
ax.text(0.003, 0.65, '100 pc', color='red')
ax.text(0.006, 0.90, '10 pc', color='orange')
ax.text(0.008, 1.02, '1 pc', color='green')

ax.add_patch(patches.Rectangle(
        (0.05, 0), 0.0912-0.05, 1.1, facecolor='blue', edgecolor='none', alpha=0.2))

cm = plt.cm.get_cmap('jet')
for i in range(1000):
    ax.add_patch(patches.Rectangle(
        (0.35 + i/1000/3, 1.02), 0.05, 0.08, facecolor=cm(i/1000), edgecolor='none'))
ax.text(0.11, 1.03, r'UV')
ax.text(0.82, 1.03, r'IR')

ax.set_ylim(0, 1.1)
ax.set_xlim(1e-4, 200)
ax.set_xscale('log')

ax.get_yaxis().set_tick_params(direction='out')
ax.get_xaxis().set_tick_params(direction='out')
ax.get_yaxis().set_tick_params(which='both', direction='out')
ax.get_xaxis().set_tick_params(which='both', direction='out')
ax.set_xlabel(r'Wavelength $\lambda$ ($\mu$m)')
ax.set_ylabel(r'Fraction of photons received')

plt.savefig(path + 'extinction.pdf', bbox_inches='tight')
#plt.show()
