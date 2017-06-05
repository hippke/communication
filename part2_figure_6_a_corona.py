import numpy
import matplotlib.pyplot as plt


data = numpy.genfromtxt('data_1961ZA.....53...81W.csv',
    skip_header=1,
    dtype=[
        ('r', 'f8'), ('d0', 'f8'), ('d10', 'f8'), ('d20', 'f8'), ('d30', 'f8'),
        ('d40', 'f8'), ('d50', 'f8'), ('d60', 'f8'), ('d70', 'f8'),
        ('d80', 'f8'), ('d90', 'f8')])


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
ax.set_xlabel(r'$b/R_{\odot}$')
ax.set_ylabel(r'Visible brightness ($L_{\odot}$)')
ax.set_xlim(1, 2)
ax.set_ylim(0.75*10**-8, 1.25*10**-6)
#ax.set_xscale('log')
ax.set_yscale('log')

y = 10**data['d0'] * 10**-12
plt.plot(data['r'], y, color='black')#, linewidth=0.5)
plt.fill_between(data['r'], y-0.05*y, y+0.05*y, color='red', alpha=0.5)
"""
y = 10**data['d10'] * 10**-12
plt.plot(data['r'], y, color='black')#, linewidth=0.5)
plt.fill_between(data['r'], y-0.05*y, y+0.05*y, color='red', alpha=0.5)

y = 10**data['d20'] * 10**-12
plt.plot(data['r'], y, color='black')#, linewidth=0.5)
plt.fill_between(data['r'], y-0.05*y, y+0.05*y, color='green', alpha=0.5)

y = 10**data['d30'] * 10**-12
plt.plot(data['r'], y, color='black')#, linewidth=0.5)
plt.fill_between(data['r'], y-0.05*y, y+0.05*y, color='c', alpha=0.5)
"""
y = 10**data['d40'] * 10**-12
plt.plot(data['r'], y, color='black')#, linewidth=0.5)
plt.fill_between(data['r'], y-0.05*y, y+0.05*y, color='blue', alpha=0.5)
"""
y = 10**data['d50'] * 10**-12
plt.plot(data['r'], y, color='black')#, linewidth=0.5)
plt.fill_between(data['r'], y-0.05*y, y+0.05*y, color='y', alpha=0.5)

y = 10**data['d60'] * 10**-12
plt.plot(data['r'], y, color='black')#, linewidth=0.5)
plt.fill_between(data['r'], y-0.05*y, y+0.05*y, color='k', alpha=0.5)

y = 10**data['d70'] * 10**-12
plt.plot(data['r'], y, color='black')#, linewidth=0.5)
plt.fill_between(data['r'], y-0.05*y, y+0.05*y, color='b', alpha=0.5)

y = 10**data['d80'] * 10**-12
plt.plot(data['r'], y, color='black')#, linewidth=0.5)
plt.fill_between(data['r'], y-0.05*y, y+0.05*y, color='green', alpha=0.5)
"""
y = 10**data['d90'] * 10**-12
plt.plot(data['r'], y, color='black')#, linewidth=0.5)
plt.fill_between(data['r'], y-0.05*y, y+0.05*y, color='green', alpha=0.5)

ax.text(1.4, 1.75*10**-7, r'$0^{\circ}$')
ax.text(1.4, 0.7*10**-7, r'$45^{\circ}$')
ax.text(1.4, 2.0*10**-8, r'$90^{\circ}$')

ay2 = ax.twiny()
ay2.get_xaxis().get_major_formatter().set_useOffset(False)
ay2.get_xaxis().get_major_formatter().set_scientific(False)
ay2.get_xaxis().set_tick_params(direction='out')
ay2.get_xaxis().set_tick_params(direction='out')
ay2.get_xaxis().set_tick_params(which='both', direction='out')
ay2.get_xaxis().set_tick_params(which='both', direction='out')
ay2.set_xlim([546, 2184])
ay2.set_xlabel(r'Heliocentric distance $z$ (au)')



plt.savefig('figure_corona.pdf', bbox_inches='tight')
# plt.show()
