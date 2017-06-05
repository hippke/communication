import numpy
import matplotlib.pyplot as plt
import matplotlib.patches as patches


path = ''
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
ax.set_ylabel(r'Achieved resolution ($D / \lambda$)')
ax.set_xlim(1E-6, 200)
ax.set_ylim(1E-7, 2)
ax.set_xscale('log')
ax.set_yscale('log')

# Diffraction Limit
plt.plot((1E-6, 200), (1, 1), color='black', linestyle='dashed')
ax.text(2 * 1E-6, 0.3, r'Diffraction limit', color='black')

# Observed limit
plt.plot((1E-6, 0.3), (3 * 1E-7, 1), color='black', linestyle='dotted')
ax.text(2 * 1E-6, 1E-3, 'Technological\nlimit 2017', color='black')

# Astrosat
plt.plot((0.00082, 0.13), (2.86615E-04, 4.54389E-02), color='blue')
plt.plot((0.13, 0.3), (4.54389E-02, 1.04859E-01), color='blue')
plt.plot((0.3, 0.53), (1.04859E-01, 1.33381E-01), color='blue')
ax.text(0.007, 0.0006, r'Astrosat', color='blue')

# Swift
plt.plot((0.000008266, 0.00008266), (1.15569E-07, 1.15569E-06), color='orange')
plt.plot((0.00008266, 0.00012398), (1.15569E-06, 5.20017E-05), color='orange')
plt.plot((0.00012398, 0.00619921), (5.20017E-05, 2.60017E-03), color='orange')
plt.plot((0.00619921, 0.17), (2.60017E-03, 2.85216E-01), color='orange')
plt.plot((0.17, 0.65), (2.85216E-01, 1.00000E+00), color='orange')
ax.text(1.5E-4, 2E-6, r'Swift', color='orange')

# Galex
plt.plot((0.1528, 0.2271), (1.78855E-02, 2.15669E-02), color='red')
ax.text(0.3, 0.008, r'Galex', color='red')

# Chandra
plt.plot((0.00012, 0.012), (5.03323E-05, 5.03323E-03), color='green')
ax.text(3 * 1E-4, 5E-5, r'Chandra', color='green')

# Spitzer
plt.plot((3.6, 4.5), (1.77643E-01, 2.22054E-01), color='purple')
plt.plot((4.5, 5.8), (2.22054E-01, 2.86203E-01), color='purple')
plt.plot((5.8, 8), (2.86203E-01, 3.94763E-01), color='purple')
plt.plot((8, 24), (3.94763E-01, 3.94763E-01), color='purple')
plt.plot((24, 70), (3.94763E-01, 5.18127E-01), color='purple')
plt.plot((70, 160), (5.18127E-01, 1.00000E+00), color='purple')
ax.text(6, 0.1, r'Spitzer', color='purple')

# Hubble
microns = (0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 2.35)
frac = (0.10066464, 0.177643482, 0.2516616, 0.3355488, 0.431419886, 0.5033232, \
    0.575226514, 0.666163059, 0.740181176, 0.814199294, 0.888217412, 0.962235529, 1)
plt.plot(microns, frac, color='brown')
ax.text(0.5, 0.30, r'HST', color='brown')

cm = plt.cm.get_cmap('jet')
for i in range(1000):
    ax.add_patch(patches.Rectangle(
        (0.35 + i/1000/3, 1*10**-7), 0.01, 2*10**-7, facecolor=cm(i/1000), edgecolor='none'))
ax.text(0.08, 1.2 * 10**-7, r'UV')
ax.text(0.90, 1.2 * 10**-7, r'IR')

plt.savefig(path + 'figure_lambda_D.pdf', bbox_inches='tight')
# plt.show()
