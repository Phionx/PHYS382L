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
# datafile = 'Analysis/20190423_CPMG_fit.csv'
# datafile = 'Analysis/20190423_T_1_fit.csv'
datafile = 'Analysis/20190423_T_2_fit.csv'
savefigs = False
with open(datafile, 'r') as csvFile:
	reader = csv.reader(csvFile)
	reader = list(reader)
csvFile.close()

concentration_set = [float(row[0]) for row in reader]
T_2_set           = [float(row[1]) for row in reader]
T_2_err_set       = [float(row[2]) for row in reader]


plt.figure(10)
plt.xlabel("Concentration (%)", fontsize=12)
plt.ylabel("$\\tau$ (s)", fontsize=12)
plt.title( "$\\tau$ vs Concentration", fontsize=14)
plt.errorbar(concentration_set, T_2_set, fmt='g')
if (savefigs):
		fig = plt.gcf()
		fig.set_size_inches((15.5, 8.5), forward=False)
		fig.savefig("figs/" + "tau vs Concentration" + ".png", dpi=500)
plt.show()