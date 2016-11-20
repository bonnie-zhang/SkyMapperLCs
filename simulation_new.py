"""

BZ 20/11/2016

SN Ia lightcurve simulation library, which:
- generates SALT2 lightcurves (using snlc) for a SN Ia, given input parameters
- simulates a SkyMapper observation of said SN
- fits the 'observed' lightcurve with snfit

Requires SALT2 to be installed (http://supernovae.in2p3.fr/salt/doku.php?id=usage). 
Set your default -S option to the location of your snlc and snfit binaries.

"""

import os
import time
import copy
import glob

from optparse import OptionParser
from collections import OrderedDict

import numpy as np
import matplotlib.pyplot as plt

## case sensitive isfile
def isfile_cs(filename):
    if not os.path.isfile(filename):
        return False
    return filename in os.listdir('.')

## distance modulus stuff

q0 = 0.2
c = 3*10**5 #km/s
H0 = 70. #km/s/Mpc

def mu(z):
    d_L = c*(z + 0.5*(1-q0)*z**2)/H0
    return 5*np.log10(d_L) + 25

filters = ['v', 'g' ,'r', 'i'][::-1] #backswards so tat legend example points are nicer colour
color_dict = {'v': 'indigo', 'g': 'forestgreen', 'r': 'r', 'i': 'gold'}

## date stamp

t = time.gmtime(time.time())
date = '%4d%02d%02d' % (t[0], t[1], t[2])
    
## set up options

parser = OptionParser()

## input parameters for SN Ia to simulate: name, date of max, extinction, stretch, colour, redshift

parser.add_option('-n', '--name', dest='name', default='test',
                  help='name of simulated SN')

parser.add_option('-d', '--daymax', dest='daymax', default=57700,
                  help='date of B-band max')

parser.add_option('-e', '--extinction', dest='MWEBV', 
                  help='MW E(B-V) of SN')

parser.add_option('-X', '--stretch', dest='X1',
                  help='SALT2 stretch X_1 of SN')

parser.add_option('-C', '--colour', dest='C', 
                  help='SALT2 colour C of SN')

parser.add_option('-z', '--redshift', dest='redshift',
                  help='redshift of SN')

## generates and renames lightcurves using snlc (don't have to do this more than once, but also it doesn't hurt)

parser.add_option('-g', '--generate', dest='generate', action="store_true", default=False,
                  help='generates snlc-simulated lightcurves')

## text file containing dates to observe SN in; see example file for format

parser.add_option('-c', '--cadenceFile', dest='cadenceFile', default = 'cadence_test.txt',
                  help='file containing dates to simulate (relative to max)')

## fit 'observed' lightcurve using SALT2 snfit

parser.add_option('-r', '--run', dest='run', action="store_true", default=False,
                  help='run snfit on simulated lightcurve')

parser.add_option('-p', '--plot', dest='plot', action="store_true", default=False,
                  help='plot lightcurves')

parser.add_option('-s', '--suffix', dest='suffix', default='',
                  help='suffix to add to observation')

## whether observed in bad seeing - this basically affects the scatter, and can probably be improved

parser.add_option('-b', '--bad_gr', dest='bad_gr', action="store_true", default=False,
                  help='use bad seeing for gr?')

parser.add_option('-i', '--bad_i', dest='bad_i', action="store_true", default=False,
                  help='use bad seeing for i?')

## bonnie has multiple salt2 versions installed so needs the following; set your default to empty string

parser.add_option('-S', '--salt', dest='salt', default='/usr/local/Cellar/salt/2.4/bin/',
                  help='location of SALT2 binaries')

(options, args) = parser.parse_args()


empty = {f: {} for f in filters}

## intrinsic scatter - can change
scatter_dict = {'v': 0.08, 'g': 0.06, 'r': 0.06, 'i': 0.1}

if options.bad_i:
    scatter_dict['i'] = 0.2 ## if observing i during bad seeing time
if options.bad_gr:
    scatter_dict['g'] = 0.12
    scatter_dict['r'] = 0.12

## name of sn
sn = options.name
z = np.float(options.redshift)

## SALT2 flux zero point - don't change this
X0 = 10**((29.69-mu(z))/2.5)

## writes simulated files: input with data, and blank lightcurve, both necessary for snlc to work

infile = open('input_' + sn + '.dat','w')
infile.write('Salt2Model\nBEGIN_OF_FITPARAMS Salt2Model\nDayMax ' + options.daymax + ' 0\nRedshift ' + str(z) +' 0 F\nColor ' + str(options.C) + ' 0\nX0 ' + str(X0) + ' 0\nX1 ' + str(options.X1) + ' 0\nEND_OF_FITPARAMS Salt2Model')
infile.close()
lcfile = open('lc_' + sn + '_.list', 'w')
lcfile.write('@Z_HELIO '+ str(z) +'\n@SN '+ sn + '\n@MWEBV ' + options.MWEBV + '\n@SNTYPE SNIa\n@SURVEY SKYMAPPER\n@INSTRUMENT SKYMAPPER\n@MAGSYS SKYMAP_AB\n#Filter:\n#Date:\n#Mag:\n#Magerr:\n#MagSys:\n')
for f in filters:
    linetowrite = 'SKYMAPPER::' + f + ' ' + ' '.join([options.daymax, '0', '0', 'SKYMAP_AB']) + '\n'
    lcfile.write(linetowrite)
lcfile.close()


## use snlc to simulate fake LCs from above input parameters; this is iterative because snlc sporadically fails

if options.generate:
    ncols = 0
    while ncols < len(filters):
        cmd_snlc = options.salt + 'snlc lc_' + sn + '_.list -p input_' + sn + '.dat'
        os.system(cmd_snlc) > 'temp'
        ## test viability of snlc outputs
        for lcfile in glob.glob('lc_Salt2Model_?.dat'):
            size = os.path.getsize(lcfile)
            print size
            if 10**3 < size < 10**4:
                cmd_mv = 'mv ' + lcfile + ' ' + lcfile[:14] + sn + '_' + lcfile[14:]
                os.system(cmd_mv)
        print glob.glob('lc_Salt2Model_' + sn + '_?.dat')
        ncols = len(glob.glob('lc_Salt2Model_' + sn + '_?.dat'))

## 'performs observation' - first set up dict of ZPs

zp_dict = {}
mag_sim, mag_sim_err, flux_sim, flux_sim_err = copy.deepcopy(empty), copy.deepcopy(empty), copy.deepcopy(empty), copy.deepcopy(empty)

lcfiles = {}

obsdates = np.genfromtxt(options.cadenceFile, names=True)
filtersObs = obsdates.dtype.names

for f in filtersObs:
    filename = 'lc_Salt2Model_' + sn + '_' + f + '.dat'
    if isfile_cs(filename): #checks if lc file exists
        lcfiles[f] = filename
for f in lcfiles.keys():
    lines = open(lcfiles[f]).readlines()
    for line in lines:
        if not (line.startswith('@') or line.startswith('#')):
            l = line.split(' ')
            try:
                if len(l) > 1 and float(l[2]) < float(l[1]):
                    date = float(l[0])
                    flux_sim[f][date] = float(l[1])
                    flux_sim_err[f][date] = float(l[2])
                    zp = float(l[3])
                    mag_sim[f][date] = -2.5*np.log10(float(l[1])) + zp
                    mag_sim_err[f][date] = 1.085736 * float(l[2])/float(l[1])
                    if f not in zp_dict.keys():
                        zp_dict[f] = zp
            except IndexError:
                pass #print '??', l

lclines = open('lc_' + sn + '_.list').readlines()


flux_obs, flux_obs_err = copy.deepcopy(empty), copy.deepcopy(empty) ## subset of flux_sim; only those 'observed'; doesn't include scatter

for f in filtersObs:
    #print f
    for p in obsdates[f]:
        if p > -99:
            d = np.float(options.daymax) + p
            i = np.argmin(abs(flux_sim[f].keys() - d))
            d_closest = flux_sim[f].keys()[i]
            flux_obs[f][d_closest] = flux_sim[f][d_closest]
            flux_obs_err[f][d_closest] = flux_sim_err[f][d_closest]

filename = 'simlc_' + sn + options.suffix + '_.list'
simfile = open(filename,'w')
simfile.write('@Z_HELIO '+ str(z) +'\n@SN '+sn + '\n@MWEBV ' + str(options.MWEBV) + '\n@SNTYPE SNIa\n@SURVEY SKYMAPPER\n@INSTRUMENT SKYMAPPER\n@MAGSYS SKYMAP_AB\n#Filter:\n#Date:\n#Mag:\n#Magerr:\n#MagSys:\n')
linestowrite = []
for f in filters:
    for d in flux_obs[f].keys():
        fl = flux_sim[f][d]
        fl_err = flux_sim_err[f][d]
        #            mag_err = np.random.normal(0,scatter_dict[f])
        mag_err = scatter_dict[f]
        mag = -2.5*np.log10(fl) + zp_dict[f] + np.random.normal(0,scatter_dict[f])
        linetowrite = 'SKYMAPPER::' + f + ' ' + ' '.join([str(d), str(mag), str(mag_err), 'SKYMAP_AB']) + '\n'
        linestowrite.append(linetowrite)
for linetowrite in linestowrite:
    simfile.write(linetowrite)
simfile.close()


if options.plot: #plots simulated from snlc and 'observed', if available
    mag_meas, mag_meas_err=copy.deepcopy(empty),copy.deepcopy(empty) ## called 'meas' to distinguish from 'obs'; even those these are the same these are read from simulated 'observed' light curves rather than sampled from snlc outputs

    try:
        lcfile = 'simlc_' + sn + '_.list'
        lines = open(lcfile).readlines()
        for line in lines:
            if not (line.startswith('@') or line.startswith('#')):
                l = line.split(' ')
                f = l[0].split(':')[2]
                date = float(l[1])
                mag_meas[f][date] = float(l[2])
                mag_meas_err[f][date] = float(l[3])
    except IOError:
        pass

    plt.figure()
    for f in filters:
        if len(mag_sim[f]) > 0:
            plt.errorbar(mag_sim[f].keys(),mag_sim[f].values(),yerr = mag_sim_err[f].values(), fmt = '+', color = color_dict[f], label = 'salt2 simulated')
        if len(mag_meas[f]) > 0:
            plt.errorbar(mag_meas[f].keys(),mag_meas[f].values(),yerr = mag_meas_err[f].values(), fmt = 'o', color = color_dict[f], label = 'measured')
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    plt.gca().invert_yaxis()
    plt.xlabel(sn + ' lightcurve')
    #plt.show()
    plt.savefig(sn + '.png')

if options.run:    
    cmd_snfit = options.salt + 'snfit -o output_' + sn + '.dat simlc_' + sn + '_.list'
    os.system(cmd_snfit) #> 'temp' ## uncomment to not see so much output
