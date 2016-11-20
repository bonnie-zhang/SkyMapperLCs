import numpy as np
import matplotlib.pyplot as plt


nSNe = 1000

fig = plt.figure()
ax = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
fig.subplots_adjust(wspace=0)
colours = ['navy', 'orange','darkturquoise', 'darkturquoise', 'forestgreen', 'forestgreen', 'crimson','crimson', 'mediumpurple', 'mediumpurple']
linestyles = ['-.','-.',':', '-', ':', '-', ':', '-',':','-']

i = 0
#for cadence in ['5ibad_v', '5ibad_vv', '10i_v', '10i_vv', '15i_v', '15i_vv', 'earlyi_vv']:
for cadence in ['4gri_5v', 'earlyiv', '3gr_10i', '3gribad','3gr_10i_v', '3gribad_v', '5gr_10i', '5gribad','5gr_10i_v', '5gribad_v']:
    fname = 'kde_data_%s_%s.txt'%(cadence, nSNe)
    try:
        kde_data = np.genfromtxt(fname, skip_header=1)
        #print cadence, kde_data[0]
        ax.plot(kde_data[:,0],kde_data[:,1], color = colours[i], linestyle = linestyles[i], label=cadence)
        ax2.plot(kde_data[:,2],kde_data[:,3], color = colours[i], linestyle = linestyles[i], label=cadence)
        i += 1
    except IOError:
        print cadence, '??'

#ax.set_xlim((-1,1))
#ax2.set_xlim((-4,4))
ax.set_xlabel(r'$C_{\rm{SALT2}}- C_{\rm{in}}$')

ax2.set_xlabel(r'$X_{1,\rm{SALT2}}- X_{1,\rm{in}}$')
ax.set_xlim((-.5,.5))
ax2.set_xlim((-2,2))
ax.axvline(0, color='k', linestyle=':', linewidth=.5)
ax2.axvline(0, color='k', linestyle=':', linewidth=.5)
ax.legend(loc=2)



fig.show()
