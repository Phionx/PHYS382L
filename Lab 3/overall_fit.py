import numpy as np
import scipy.optimize as opt
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import math
import csv
import scipy.special as special


# For displaying values
from decimal import Decimal

#IMPORT DATA
#-----------------------------------------------------------------------------------------------------------------------------------------
# datafile = 'Analysis/20190424_T_1_fit.csv'
# ylabel = "$T_1$ (s)"
# title  = "$T_1$ vs Concentration"

# datafile = 'Analysis/20190424_CPMG_fit.csv'
# ylabel = "$T_2$ (s)"
# title  = "$T_2$ vs Concentration"

datafile = 'Analysis/20190424_CPMG_fit_T_2_star.csv'
ylabel = "$T_2^*$ (s)"
title  = "$T_2^*$ vs Concentration"

savefigs = False
with open(datafile, 'r') as csvFile:
	reader = csv.reader(csvFile)
	reader = list(reader)
csvFile.close()

concentration_set = np.array([float(row[0]) for row in reader])
T_2_set           = np.array([float(row[1]) for row in reader])
T_2_err_set       = np.array([float(row[2]) for row in reader])


print(concentration_set)

plt.figure(10)
# plt.ylim(0, max(T_2_set))
plt.xlabel("Concentration (%)", fontsize=12)
plt.ylabel(ylabel, fontsize=12)
plt.title( title, fontsize=14)
# plt.errorbar(concentration_set, T_2_set, T_2_err_set, fmt='g')

plt.plot(concentration_set, T_2_set, '--ko', markersize=3, color='r')#color='#3F7F4C')
plt.fill_between(concentration_set, T_2_set-T_2_err_set, T_2_set+T_2_err_set, alpha=.8, edgecolor='#3F7F4C', facecolor='#7EFF99', linewidth=0)

if (savefigs):
		fig = plt.gcf()
		fig.savefig("figs/" + title + ".png", dpi=500)
else:
	plt.show()