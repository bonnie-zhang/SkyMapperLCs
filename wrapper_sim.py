"""

Wrapper for simulation_new.py. BZ wrote this to assess the recovered SALT2 X1 and C for a number of different observing cadences, for the same 1000 SNe. This is here as an example: play with and change it as you see fit, and send me questions about parts you don't understand!

"""

# pylint: disable=E1101

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

from scipy.optimize import curve_fit

from KDE_utils import solve_kde

seed = 999 

np.random.seed(seed)

## edit below for simulation:

nSNe = 5 ## change this for more SNe

## generate random phase of first detection, and time to observe over
tDet = np.random.randint(-15, -2, size=nSNe)
tObs = np.random.randint(25, 65, size=nSNe)

## generate random SN properties, redshifts, extinction
z_in = np.random.uniform(0.01, 0.08, size=nSNe)
X1_in = 3*np.random.randn(nSNe)
C_in = 0.3*np.random.randn(nSNe)
MWEBV_in = 0.008*np.random.randn(nSNe) + 0.01
daymax_in = np.random.uniform(57650, 57700, size=nSNe) #doesn't really matter

## toggle for observing gr in bad seeing time, and for including v observations around max

#bad_gr = True
"""

cad_gr = 5

bad_i = False
cad_i = 10

Nv = 0 ## n of v observations observations

"""


## for plotting KDEs

def gauss_fn(x,mu,sigma,A): #gaussian for fitting to the probability distributions
    return A*np.exp(-(x-mu)**2/(2*sigma**2))


def bandwidth(data):
    return 1.06 * np.std(data) / np.power(len(data),1 / 5)


for cad_gr in [5]:#,3]:
    for Nv in [0]:#,1]:
        for cad_i in [cad_gr]:#, 10]:
            if cad_i == 10:
                bad_i = False
            else:
                bad_i = True
                
            ## name of cadence, writes to file below
            if cad_gr == cad_i:
                cadence = '%sgribad'%cad_gr + '_v'*Nv
            else:
                cadence = '%sgr_%si'%(cad_gr, cad_i) + '_v'*Nv

            print cadence
            
            ## wrapper for simulation for each SN

            if Nv > 0:
                filtersObs = ['g','r','i','v']
            else:
                filtersObs = ['g','r','i']
            for i in range(nSNe):
                ## record cadence as array
                nObs = len(np.arange(tDet[i],tDet[i]+tObs[i],cad_gr))
                arr = -99*np.ones((nObs,len(filtersObs)))
                arr[:,0] = np.arange(tDet[i],tDet[i]+tObs[i],cad_gr)
                arr[:,1] = np.arange(tDet[i],tDet[i]+tObs[i],cad_gr)
                l = len(np.arange(tDet[i],tDet[i]+tObs[i],cad_i))
                arr[:l,2] = np.arange(tDet[i],tDet[i]+tObs[i],cad_i)
                if Nv == 1:
                    arr[:1,3] = [2]
                name = '%s_%s'%(cadence,i)
                fhandle = open('cadence_' + name + '.txt','w+')
                fhandle.write(' '.join(filtersObs) + '\n')
                for k in range(nObs):
                    fhandle.write(' '.join([str(int(x)) for x in arr[k]]) + '\n') 
                fhandle.close()
                ## only use flags -g, -r if necessary: if outputs generated, can comment out these flags to save time (having seeded a random number allows this)
                if bad_i:
                    cmd = 'python simulation_new.py -n %s -d %s -e %s -X %s -C %s -z %s -c cadence_%s.txt -p -r -g -b -i'%(name, daymax_in[i], MWEBV_in[i], X1_in[i], C_in[i], z_in[i], name)
                else:
                    cmd = 'python simulation_new.py -n %s -d %s -e %s -X %s -C %s -z %s -c cadence_%s.txt -p -r -g -b'%(name, daymax_in[i], MWEBV_in[i], X1_in[i], C_in[i], z_in[i], name)
                #print cmd
                os.system(cmd)
                
            ## document failed SALT2 fits, and output quantities
            fail = 0

            mB_out = []
            mB_out_err = []
            X1_out = []
            X1_out_err = []
            C_out = []
            C_out_err = []

            for i in range(nSNe):
                fname = 'output_%s_%s.dat'%(cadence,i)
                try:
                    lines = open(fname).readlines()
                    for line in lines:
                        if line.startswith('Rest'):
                            mB_out.append(float(line.split()[1]))
                            mB_out_err.append(float(line.split()[2]))
                            break
                    for line in lines:
                        if line.startswith('Color'):
                            C_out.append(float(line.split()[1]))
                            C_out_err.append(float(line.split()[2]))
                            break
                    for line in lines:
                        if line.startswith('X1'):
                            X1_out.append(float(line.split()[1]))
                            X1_out_err.append(float(line.split()[2]))
                            break
                except IOError:
                    print '? %s not found'%name
                    fail += 1
                    mB_out.append(0)
                    mB_out_err.append(0)
                    C_out.append(0)
                    C_out_err.append(0)
                    X1_out.append(0)
                    X1_out_err.append(0)


            fig = plt.figure()

            ## not too much point plotting m_B, so commented out below
            #gs = gridspec.GridSpec(3, 2, width_ratios=[3,1])
            gs = gridspec.GridSpec(2, 2, width_ratios=[3,1])
            gs.update(wspace=0, hspace=0)
            
            ## colour
            
            ax_C = plt.subplot(gs[0,0])
            ax_C.plot(range(nSNe), C_in, linestyle='None', marker='^', color = 'm', label='input colour')
            ax_C.errorbar(range(nSNe), C_out, yerr = C_out_err,linestyle='None', fmt = '+', mfc='navy', ecolor='navy', label='simulated with cadence %s'%cadence)
            ax_C.set_ylabel(r'$C$')
            ax_C.set_xlim((-0.5,nSNe - 0.5))
            
            w = np.ones(nSNe)/nSNe
            ax_C_hist = plt.subplot(gs[0,1], sharey = ax_C)
            ax_C_hist.hist(C_in, histtype='step', color='m', bins = np.arange(-1,1,0.05), normed = 0, weights = w, orientation='horizontal')
            ax_C_hist.hist(C_out, histtype='stepfilled', color='navy', alpha = 0.5, bins = np.arange(-1,1,0.05), normed = 0, weights = w, orientation='horizontal')
            ax_C.set_ylim((C_in.min(), C_in.max()))
            plt.setp(ax_C_hist.get_yticklabels(), visible=False)
            plt.setp(ax_C_hist.get_xticklabels(), visible=False)
            
            ## stretch
            
            ax_X1 = plt.subplot(gs[1,0], sharex = ax_C)
            ax_X1.plot(range(nSNe), X1_in, linestyle='None', marker='^', color = 'm', label='input colour')
            ax_X1.errorbar(range(nSNe), X1_out, yerr = X1_out_err,linestyle='None', fmt = '+', mfc='navy', ecolor='navy', label='simulated with cadence %s'%cadence)
            ax_X1.set_ylabel(r'$X_1$')
            ax_X1.set_xlim((-0.5,nSNe - 0.5))
            
            w = np.ones(nSNe)/nSNe
            ax_X1_hist = plt.subplot(gs[1,1], sharey = ax_X1)
            ax_X1_hist.hist(X1_in, histtype='step', color='m', bins = np.arange(-7,7,0.5), normed = 0, weights = w, orientation='horizontal')
            ax_X1_hist.hist(X1_out, histtype='stepfilled', color='navy', alpha = 0.5, bins = np.arange(-7,7,0.5), normed = 0, weights = w, orientation='horizontal')
            ax_X1.set_ylim((X1_in.min(), X1_in.max()))
            plt.setp(ax_X1_hist.get_yticklabels(), visible=False)
            plt.setp(ax_X1_hist.get_xticklabels(), visible=False)


            ax_X1.set_xlabel('simulation number')
            ax_C.set_title('simulated stretch/colour with cadence %s'%cadence)
   

            fig.show()
            
            ## save figure
            plt.savefig('results_%s_%s.png'%(cadence, nSNe))


            """## difference, with KDE
            
            fig2 = plt.figure()
            
            C_diff = C_out - C_in
            X1_diff = X1_out - X1_in

            gs = gridspec.GridSpec(1, 2)
            gs.update(wspace=0)

            ## colour

            ax = plt.subplot(gs[0,0])

            ax.hist(C_diff, histtype='step', color='m', bins = np.arange(C_diff.min(), C_diff.max(), 0.05), normed=0, weights=np.ones_like(C_diff)/len(C_diff))

            xlist = np.linspace(min(C_diff) - 2,max(C_diff) + 2,300)
            kde_array = solve_kde(xlist,C_diff,np.array(C_out_err),bandwidth(C_diff))
            kde_vector_data = np.sum(np.sum(kde_array,axis=2),axis=1)
            
            ax.plot(xlist,kde_vector_data/nSNe*0.05,color='red',label='KDE')
            ax.set_xlabel(r'$C_{\rm{SALT2}}- C_{\rm{in}}$')
            ax.set_ylim((0, kde_vector_data.max()/nSNe*0.05))
            ax.set_xlim((-1,1))

            ax2 = plt.subplot(gs[0,1])
            ax2.hist(X1_diff, histtype='step', color='m', bins = np.arange(X1_diff.min(), X1_diff.max(), 0.5), normed=0, weights=np.ones_like(X1_diff)/len(X1_diff))

            xlist2 = np.linspace(min(X1_diff) - 2,max(X1_diff) + 2,300)
            kde_array2 = solve_kde(xlist2,X1_diff,np.array(X1_out_err),bandwidth(X1_diff))
            kde_vector_data2 = np.sum(np.sum(kde_array2,axis=2),axis=1)
            ax2.plot(xlist2,kde_vector_data2/nSNe*0.5,color='red',label='KDE')
            ax2.set_xlabel(r'$X_{1,\rm{SALT2}}- X_{1,\rm{in}}$')
            ax2.set_ylim((0, kde_vector_data2.max()/nSNe*0.5))
            ax2.set_xlim((-4,4))
            #plt.title('differences in stretch/colour from simulation with cadence %s'%cadence)
            
            plt.show()
            plt.savefig('results_diff_%s_%s.png'%(cadence, nSNe))

            ## write
            
            kde_data = np.column_stack((xlist, kde_vector_data/nSNe*.05, xlist2, kde_vector_data2/nSNe*.5))
                        
            np.savetxt('kde_data_%s_%s.txt'%(cadence, nSNe), kde_data, header = 'diff_C, KDE_C, diff_X1, KDE_X1')

            
            ## fit gaussians to 
            popt_C, pcov_C = curve_fit(gauss_fn, xlist, kde_vector_data, p0 = [0, 1, 1])
            popt_X1, pcov_X1 = curve_fit(gauss_fn, xlist2, kde_vector_data2, p0 = [0, 1, 1])

            print 'color KDE: mean %s, stddev %s'%(popt_C[0], popt_C[1])
            print 'stretch KDE: mean %s, stddev %s'%(popt_X1[0], popt_X1[1])


            ## save figures of merit
            resfhandle = open('results_%s_%s.txt'%(cadence, nSNe),'w')
            resfhandle.write(cadence + '\n')
            resfhandle.write('failed: %s \n'%fail)
    
            chisqC = 0
            chisqX1 = 0
            
            for i in range(nSNe):
                if C_out_err[i] != 0:
                    chisqC += (C_in[i] - C_out[i])**2/C_out_err[i]**2
                if X1_out_err[i] != 0:
                    chisqX1 += (X1_in[i] - X1_out[i])**2/X1_out_err[i]**2
            
            redchisqC = chisqC/(nSNe-fail)
            redchisqX1 = chisqX1/(nSNe-fail)
            
            print redchisqC, redchisqX1, fail

            resfhandle.write('colour chi^2/nSNe: %s \n'%redchisqC)
            resfhandle.write('stretch chi^2/nSNe: %s \n'%redchisqX1)
            
            resfhandle.write('color KDE: mean %s, stddev %s \n'%(popt_C[0], popt_C[1]))
            resfhandle.write('stretch KDE: mean %s, stddev %s \n'%(popt_X1[0], popt_X1[1]))
            resfhandle.close()"""
