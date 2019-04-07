import numpy as np
import scipy.optimize as opt
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import math
import csv


# For displaying values
from decimal import Decimal

#DATA PARSING
#-----------------------------------------------------------------------------------------------------------------------------------------
#T,E,M,C,X  - Temperature, Energy, Magnetism, Specific Heat, Susceptibility
EMfile = 'Data/NDep/4_05_N50_1.20T3.00_EM_v0.csv'
with open(EMfile, 'r') as csvEMFile:
	readerEM = csv.reader(csvEMFile)
	readerEM   = list(readerEM)
csvEMFile.close()
readerEM  = readerEM[4:]
# first_col = [row[0] for row in readerEM]
# max_ind   = first_col.index('')
# readerEM  = readerEM[:max_ind]

T_E    = []
E_mean = []
E_std  = []

T_M    = []
M_mean = []
M_std  = []

T_C    = []
C_mean = []
C_std  = []

T_X    = []
X_mean = []
X_std  = []

for row in readerEM:
	if row[2] != "":
		T_E.append(row[0])
		E_mean.append(row[2])
		E_std.append(row[3])

	if row[4] != "":
		T_M.append(row[0])
		M_mean.append(row[4])
		M_std.append(row[5])

	if row[6] != "":
		T_C.append(row[0])
		C_mean.append(row[6])
		C_std.append(row[7])

	if row[8] != "":
		T_X.append(row[0])
		X_mean.append(row[8])
		X_std.append(row[9])

T_E    = np.array(T_E).astype(np.float)
E_mean = np.array(E_mean).astype(np.float)
E_std  = np.array(E_std).astype(np.float)

T_M    = np.array(T_M).astype(np.float)
M_mean = np.array(M_mean).astype(np.float)
M_mean = np.abs(M_mean)
M_std  = np.array(M_std).astype(np.float)

T_C    = np.array(T_C).astype(np.float)
C_mean = np.array(C_mean).astype(np.float)
C_std  = np.array(C_std).astype(np.float)


T_X    = np.array(T_X).astype(np.float)
X_mean = np.array(X_mean).astype(np.float)
X_std  = np.array(X_std).astype(np.float)

print("Done Processing EM Data\n")



#Spin Correlations
SCfile = 'Data/APR:3:2019/CLEANED_3_30_N100_1.20T3.00_SC_v0.csv'
with open(SCfile, 'r') as csvSCFile:
	readerSC = csv.reader(csvSCFile)
	readerSC = list(readerSC)
csvSCFile.close()
d_max     = int(readerSC[1][0])/2-1 #N value/2
readerSC  = readerSC[4:]
# first_col = [row[0] for row in readerSC]
# readerSC  = readerSC[:max_ind]

#Data Parsing
T_SC    = []
SC_mean = []
SC_std  = []

for row in readerSC:
	if row[2] != "":
		T_SC.append(row[0])
		SC_mean.append(np.array(row[2:d_max+2]).astype(np.float))
		SC_std.append(np.array(row[d_max+2:2*d_max+2]).astype(np.float))

T_SC    = np.array(T_SC).astype(np.float)
SC_mean = np.array(SC_mean).astype(np.float)
SC_std  = np.array(SC_std).astype(np.float)
print("Done Processing SC Data\n")
#-----------------------------------------------------------------------------------------------------------------------------------------

#FUNCTIONS
#-----------------------------------------------------------------------------------------------------------------------------------------
#Chi Square Fit

def chi_square(function, xdata, ydata, yerror, *parameter_vals):
	ytheory         = [0 for x in range(len(xdata))]
	for i in range(len(xdata)):
		ytheory[i]  = function(xdata[i], *parameter_vals)
	
	total_residuals = 0
	length          = len(yerror)
	for i in range(length):
		total_residuals = total_residuals + ((ytheory[i] - ydata[i])**2)/(yerror[i]**2)
	if (len(parameter_vals) == 0):
		return total_residuals
	
	total_residuals = total_residuals/len(parameter_vals)
	return total_residuals



def optimal_fit(function, xdata, ydata, yerror, parameter_estimates, parameter_ranges):
	popt, pcov = opt.curve_fit(function, xdata, ydata, sigma =yerror, p0=parameter_estimates,  bounds=parameter_ranges)
	perr = np.sqrt(np.diag(pcov))
	return (popt, perr)

#FITS: #t = (T - T_c)/T_c
#Energy Full Fit
#-------------------------------------------------------------------------
def kappa_calc(x, *parameters):
	T     = x
	return 2.0*np.sinh(2.0/T)/((np.cosh(2.0/T))**2.0) 

def K_one_calc(x, *parameters):
	kappa = x
	return integrate.quad(lambda phi: (1-(kappa**2.0)*((np.sin(phi))**2))**(-1.0/2.0), 0.0, np.pi/2.0)[0]

def energy_func(x, *parameters):
	# scale = parameters[0]
	# add   = parameters[1]

	#number
	if (np.isscalar(x)):
		T = x
		kappa  = kappa_calc(T)
		k_one  = K_one_calc(kappa)
		answer =  -2.0*np.tanh(1.0/T) - (((np.sinh(2.0/T))**2.0 - 1.0)/((np.sinh(2.0/T))*(np.cosh(2.0/T))))*((2.0/np.pi)*k_one - 1.0)
		return answer

	#Array
	T_set = x
	ans   = [0 for T in T_set]
	data_points = len(ans)
	for i in range(data_points):
		T = T_set[i]
		kappa  = kappa_calc(T)
		k_one  = K_one_calc(kappa)
		answer =  -2.0*np.tanh(1.0/T) - (((np.sinh(2.0/T))**2.0 - 1.0)/((np.sinh(2.0/T))*(np.cosh(2.0/T))))*((2.0/np.pi)*k_one - 1.0)
		ans[i] = answer
	ans = np.array(ans)
	return ans

#-------------------------------------------------------------------------

def mag_func(x, *parameters):
	T     = x
	ans = (1-(np.sinh(2.0/T))**(-4.0))**(1.0/8.0)
	if math.isnan(ans):
		return 0
	return ans

#-------------------------------------------------------------------------

#Theory  Value: 0
def alpha_func(x, *parameters):
	T     = x
	T_c   = parameters[0]
	alpha = parameters[1]
	constant = parameters[2]

	y = constant*np.abs((T - T_c)/T_c)**(-1.0*alpha)
	return y

#-------------------------------------------------------------------------

#Magnetization: |M| \propto |t|^{\beta}
#Theory  Value: 1/8
def beta_func(x, *parameters):
	T     = x
	T_c   = parameters[0]
	beta = parameters[1]
	constant = parameters[2]

	y = constant*np.abs((T - T_c)/T_c)**(1.0*beta)
	return y

#-------------------------------------------------------------------------

#Susceptibility: \chi \propto |t|^{-\gamma}
#Theory  Value: 7/4
def gamma_func(x, *parameters):
	T     = x
	T_c   = parameters[0]
	gamma = parameters[1]
	constant = parameters[2]

	y = constant*np.abs((T - T_c)/T_c)**(-1.0*gamma)
	return y

#-------------------------------------------------------------------------

#\xi \propto |t|^{-\nu}
#Theory  Value: 1
def nu_func(x, *parameters):
	T     = x
	T_c   = parameters[0]
	nu    = parameters[1]
	constant = parameters[2]

	y = constant*np.abs((T - T_c)/T_c)**(-1.0*nu)
	return y

#-------------------------------------------------------------------------

#R(x) = |x|^{-\eta}
#Theory  Value: 1/4
def eta_func(x, *parameters):
	T_c   = parameters[0]
	eta   = parameters[1]
	constant = parameters[2]

	y = constant*np.abs(x)**(-1.0*eta)
	return y

#-------------------------------------------------------------------------
savefigs = True
#Plotting
def plot_fit(xdata, ydata, yerror, xtheory, ytheory, params, params_err, params_names, fig_num, **graph_labels):
	x_label = graph_labels['x_label']
	y_label = graph_labels['y_label']
	title   = graph_labels['title']
	fit_eq  = graph_labels['fit_eq']

	xdata   = np.array(xdata)
	ydata   = np.array(ydata)
	yerror  = np.array(yerror)
	xtheory = np.array(xtheory)
	ytheory = np.array(ytheory)
	
	#Experiment
	plt.figure(fig_num)
	plt.plot(xdata, ydata, 'ko', markersize=1, color='#3F7F4C', label="Simulation Data")
	plt.fill_between(xdata, ydata-yerror, ydata+yerror, alpha=.8, edgecolor='#3F7F4C', facecolor='#7EFF99', linewidth=0)

	#Best Theory
	theory_label = "Theory"

	
	theory_label += " " + fit_eq

	for i in range(len(params)):
		theory_label += "\n"
		theory_label += params_names[i] + ": "
		theory_label += ("{:.2E}".format(Decimal(params[i])))
		if str(params_err[i]) != "":
			theory_label += " $\\pm$ " + (("{:.2E}".format(Decimal(params_err[i]))))

	plt.plot(xtheory, ytheory, 'k-', color='#CC4F1B', label=theory_label)
	plt.legend(loc='upper right')
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	plt.title(title)
	if (savefigs):
		fig = plt.gcf()
		fig.set_size_inches((14, 9.5), forward=False)
		fig.savefig("figs/" + title+".png", dpi=500)

	


def plot_residuals(xdata, ydata, yerror, xtheory, ytheory, fig_num, **graph_labels):
	x_label = graph_labels['x_label']
	y_label = graph_labels['y_label']
	title   = graph_labels['title']

	xdata   = np.array(xdata)
	ydata   = np.array(ydata)
	yerror  = np.array(yerror)
	xtheory = np.array(xtheory)
	ytheory = np.array(ytheory)
	data_points = len(xdata)

	#RESIDUALS
	print("xdata len: " + str(len(xdata)) + " " + "ydata len: " + str(len(ydata)))
	# print("residuals(data): " + str(ydata))
	# print("residuals(theory): " + str(ytheory))
	yres    = [(ydata[i] - ytheory[i])/yerror[i] for i in range(data_points)]
	xres    = xdata
	plt.figure(fig_num)
	plt.plot(xres, yres, 'k-', color='#3F7F4C')
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	plt.title(title + " (Residuals)")

	if (savefigs):
		fig = plt.gcf()
		fig.set_size_inches((14, 9.5), forward=False)
		fig.savefig("figs/" + title + " (Residuals).png", dpi=400)

	#INSET
	plt.figure(fig_num-1)
	try:
		if(graph_labels['log_scale_x'] and graph_labels['log_scale_y']):
				plt.loglog(xdata, ydata, 'ko', markersize=1, color='#3F7F4C', label="Simulation Data")
				# plt.fill_between(xdata, ydata-yerror, ydata+yerror, alpha=.8, edgecolor='#3F7F4C', facecolor='#7EFF99', linewidth=0)
				plt.loglog(xtheory, ytheory, 'k-', color='#CC4F1B')
	except:
		plt.plot(xdata, ydata, 'ko', markersize=1, color='#3F7F4C', label="Simulation Data")
		plt.fill_between(xdata, ydata-yerror, ydata+yerror, alpha=.8, edgecolor='#3F7F4C', facecolor='#7EFF99', linewidth=0)
		plt.plot(xtheory, ytheory, 'k-', color='#CC4F1B')
	
	# plt.legend(loc='upper right')
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	plt.title(title + " (Inset)")
	

	if (savefigs):
		fig = plt.gcf()
		fig.set_size_inches((14, 9.5), forward=False)
		fig.savefig("figs/" + title + " (Inset).png", dpi=500)

	return

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
#-----------------------------------------------------------------------------------------------------------------------------------------


#Fitting
#-----------------------------------------------------------------------------------------------------------------------------------------
#Estimates




# Estimate_T_c_nu      = 2.2
# Estimate_nu          = 1.0

# Estimate_eta         = 1.0/4



#Energy Full Fit
#-----------------------------------------------------------------------------------------------------------------------------------------
print("Energy Full Fit: BEGIN\n")
scale_constant            = 1.0
add_constant              = 1.0
parameter_estimates       = []

var                       = [[],[]]
parameter_bound_ranges    = ([parameter_estimates[i] + var[0][i] for i in range(len(var[0]))], [parameter_estimates[i] + var[1][i] for i in range(len(var[1]))])
xdata_data                = T_E
ydata_data                = E_mean
yerror_data               = E_std
function                  = energy_func
fig_num                   = 1

xdata_fit                 = xdata_data
ydata_fit                 = ydata_data
yerror_fit                = yerror_data



Optimal_params        = []
Optimal_params_err    = []
Optimal_params_names  = []



xfit = xdata_fit
yfit = function(xfit)

title    = "Analytic Fit of Energy vs. Temperature (N = 50)"
y_label  = "Energy (J)"
x_label  = "Temperature (J/$k_B$)"
fit_eq   = "$\\hat{E}(T) = A*E(T) + B$"

plot_fit(xdata_data, ydata_data, yerror_data, xfit, yfit, Optimal_params, Optimal_params_err, Optimal_params_names, fig_num, title=title, x_label=x_label, y_label=y_label, fit_eq=fit_eq)

fig_num += 2
# plot_residuals(xdata_fit, ydata_fit, yerror_fit, xfit, yfit, fig_num, title=title, x_label=x_label, y_label=y_label)

print("Energy Full Fit: END\n")
plt.show()
#-----------------------------------------------------------------------------------------------------------------------------------------
