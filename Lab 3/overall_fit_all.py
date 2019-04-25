import numpy as np
import scipy.optimize as opt
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import math
import csv
import scipy.special as special
from random import randint


# For displaying values
from decimal import Decimal

#IMPORT DATA
#-----------------------------------------------------------------------------------------------------------------------------------------
datafile_T_1      = 'Analysis/20190425_T_1_fit_visc.csv'
datafile_T_2_star = 'Analysis/20190425_CPMG_fit_T_2_star_visc.csv'
datafile_T_2      = 'Analysis/20190425_CPMG_fit_visc.csv'

datafiles = [datafile_T_1, datafile_T_2_star, datafile_T_2]
labels    = ["$T_1$", "$T_2^*$", "$T_2$", "$\\frac{1}{1/(2T_1) + 1/T_2^*}$"]
ylabel = "$\\tau$ (s)"


savefigs = True
readers = [None for x in range(3)]

for i in range(3):
	with open(datafiles[i], 'r') as csvFile:
		readers[i] = csv.reader(csvFile)
		readers[i] = list(readers[i])
	csvFile.close()


concentration_set = np.array([50.0, 60.0, 70.0, 80.0, 90.0])
visc_set          = np.array([.00771, .0146, .0315, .0803, .260])

T_set     = [None for x in range(4)]
T_err_set = [None for x in range(4)]

for i in range(3):
	T_set[i]           = np.array([float(row[1]) for row in readers[i]])
	T_err_set[i]       = np.array([float(row[2]) for row in readers[i]])



def combo(T_1, T_2_star):
	return 1.0/(1.0/(2.0*T_1) + 1.0/T_2_star)

def combo_err(T_1, T_1_err, T_2_star, T_2_star_err):
	ans  =  (T_1_err**(2.0))*(1.0/(2.0*T_1**2.0)*combo(T_1, T_2_star)**2.0)**2.0
	ans  += (T_2_star_err**(2.0))*(1.0/(T_2_star**2.0)*combo(T_1, T_2_star)**2.0)**2.0
	ans  = np.sqrt(ans)
	return ans

T_set[3]     = np.array([combo(T_set[0][i], T_set[1][i]) for i in range(len(T_set[0]))])
T_err_set[3] = np.array([combo_err(T_set[0][i], T_err_set[0][i], T_set[1][i], T_err_set[1][i]) for i in range(len(T_set[0]))])

# T_set[0]     = 2.0*T_set[0]
# T_err_set[0] = 2.0*T_err_set[0]

#Viscosity Graph
#-----------------------------------------------------------------------------------------------------------------------------------------
plt.figure(1)
title  = "$\\tau$ vs Viscosity"
plt.xlabel("Viscosity ($Ns/m^2$)", fontsize=12)
plt.ylabel(ylabel, fontsize=12)
plt.title( title, fontsize=14)


colors    = []
for i in range(8):
    colors.append('#%06X' % randint(0, 0xFFFFFF))

i=0
# plt.errorbar(visc_set, T_set[i], T_err_set[i], fmt='--ro', markersize=3,  label=labels[i])
plt.plot(visc_set, T_set[i], '--mo', markersize=3, label=labels[i])#color='#3F7F4C', color='r')
plt.fill_between(visc_set, T_set[i]-T_err_set[i], T_set[i]+T_err_set[i], alpha=.8,  edgecolor='#3F7F4C', facecolor='#7EFF99', linewidth=0) #edgecolor='#3F7F4C', facecolor='#7EFF99',

i=2
# plt.errorbar(visc_set, T_set[i], T_err_set[i], fmt='--bo', markersize=3,  label=labels[i])
plt.plot(visc_set, T_set[i], '--ko', markersize=3, label=labels[i])#color='#3F7F4C', color='r')
plt.fill_between(visc_set, T_set[i]-T_err_set[i], T_set[i]+T_err_set[i], alpha=.8,  edgecolor='r', facecolor='r', linewidth=0) #edgecolor='#3F7F4C', facecolor='#7EFF99',

# for i in [0,2]:
	# plt.plot(visc_set, T_set[i], '--ko', markersize=3, color=colors[2*i], label=labels[i])#color='#3F7F4C', color='r')
	# plt.fill_between(visc_set, T_set[i]-T_err_set[i], T_set[i]+T_err_set[i], alpha=.8,  edgecolor=colors[2*i+1], facecolor=colors[2*i+1], linewidth=0) #edgecolor='#3F7F4C', facecolor='#7EFF99',

plt.legend(loc=1, prop={'size': 8.5})

if (savefigs):
		fig = plt.gcf()
		fig.savefig("figs/" + title + ".png", dpi=500)
else:
	plt.show()
#-----------------------------------------------------------------------------------------------------------------------------------------


#Concentration Graph
#-----------------------------------------------------------------------------------------------------------------------------------------
plt.figure(2)
title  = "$\\tau$ vs Concentration"
plt.xlabel("Concentration (%)", fontsize=12)
plt.ylabel(ylabel, fontsize=12)
plt.title( title, fontsize=14)


colors    = []
for i in range(8):
    colors.append('#%06X' % randint(0, 0xFFFFFF))

i=0
# plt.errorbar(concentration_set, T_set[i], T_err_set[i], fmt='--ro', markersize=3,  label=labels[i])
plt.plot(concentration_set, T_set[i], '--mo', markersize=3, label=labels[i])#color='#3F7F4C', color='r')
plt.fill_between(concentration_set, T_set[i]-T_err_set[i], T_set[i]+T_err_set[i], alpha=.8,  edgecolor='#3F7F4C', facecolor='#7EFF99', linewidth=0) #edgecolor='#3F7F4C', facecolor='#7EFF99',

i=2
# plt.errorbar(concentration_set, T_set[i], T_err_set[i], fmt='--bo', markersize=3,  label=labels[i])
plt.plot(concentration_set, T_set[i], '--ko', markersize=3, label=labels[i])#color='#3F7F4C', color='r')
plt.fill_between(concentration_set, T_set[i]-T_err_set[i], T_set[i]+T_err_set[i], alpha=.8,  edgecolor='r', facecolor='r', linewidth=0) #edgecolor='#3F7F4C', facecolor='#7EFF99',

# for i in [0,2]:
	# plt.plot(visc_set, T_set[i], '--ko', markersize=3, color=colors[2*i], label=labels[i])#color='#3F7F4C', color='r')
	# plt.fill_between(visc_set, T_set[i]-T_err_set[i], T_set[i]+T_err_set[i], alpha=.8,  edgecolor=colors[2*i+1], facecolor=colors[2*i+1], linewidth=0) #edgecolor='#3F7F4C', facecolor='#7EFF99',

plt.legend(loc=1, prop={'size': 8.5})

if (savefigs):
		fig = plt.gcf()
		fig.savefig("figs/" + title + ".png", dpi=500)
else:
	plt.show()