import numpy as np
import scipy.optimize as opt
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import math
import csv
import scipy.special as special


# For displaying values
from decimal import Decimal

#DATA PARSING
#-----------------------------------------------------------------------------------------------------------------------------------------
#T,E,M,C,X  - Temperature, Energy, Magnetism, Specific Heat, Susceptibility
datafiles      = ['Data/20190418/CPMG_s5.csv','Data/20190418/CPMG_s6.csv', 'Data/20190418/CPMG_s7.csv', 'Data/20190418/CPMG_s8.csv', 'Data/20190418/CPMG_s9.csv']
concentrations = [.50, .60, .70, .80, .90]
concentrations = [x*100 for x in concentrations]
num_files   = len(datafiles)

times       = []
sigma_Y     = []
times_err   = []
sigma_Y_err = []
for i in range(num_files):
	datafile   = datafiles[i]
	with open(datafile, 'r') as csvFile:
		reader = csv.reader(csvFile)
		reader = list(reader)
		reader = reader[1:]
	csvFile.close()

	concentration = concentrations[i] #concentration

	times_i       = [row[0] for row in reader]
	sigma_Y_i     = [row[1] for row in reader]
	times_err_i   = [row[2] for row in reader]
	sigma_Y_err_i = [row[3] for row in reader]

	times_i       = np.array(times_i).astype(np.float)
	times_err_i   = np.array(times_err_i).astype(np.float)

	sigma_Y_i     = np.array(sigma_Y_i).astype(np.float)
	sigma_Y_err_i = np.array(sigma_Y_err_i).astype(np.float)
	
	times.append(times_i)
	times_err.append(times_err_i)

	sigma_Y.append(sigma_Y_i)
	sigma_Y_err.append(sigma_Y_err_i)

times       = np.array(times)
times_err   = np.array(times_err)
sigma_Y     = np.array(sigma_Y)
sigma_Y_err = np.array(sigma_Y_err) 


print("Done Processing Data\n")
#-----------------------------------------------------------------------------------------------------------------------------------------

#CHI SQUARE
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
	# popt, pcov = opt.curve_fit(function, xdata, ydata, sigma =yerror, p0=parameter_estimates,  bounds=parameter_ranges)
	# popt, pcov = opt.curve_fit(function, xdata, ydata, p0=parameter_estimates, sigma=yerror, maxfev=10**7)
	popt, pcov = opt.curve_fit(function, xdata, ydata, p0=parameter_estimates, bounds=parameter_ranges, sigma=yerror, maxfev=10**7)
	perr = np.sqrt(np.diag(pcov))
	return (popt, perr)

#FIT EQUATIONS
#-----------------------------------------------------------------------------------------------------------------------------------------
#PI PULSE (T_2^*)

def CPMG_func(x, *parameters):
	#tau, delta, t_0, A, phi, D
	t      = x
	tau    = parameters[0] 
	t_0    = parameters[1]
	A      = parameters[2]
	D      = parameters[3]

	ans    = -1.0*A*np.e**(-1.0*(t-t_0)/tau) + D
	return ans
	
def pi_2_pulse_err_func(x, xerr, *parameters): #the same as pi_pulse_err_func cause t = t_0
	t      = x
	tau    = parameters[0] 
	t_0    = parameters[1]
	A      = parameters[2]
	D      = parameters[3]

	t_err      = xerr
	tau_err    = parameters[4]
	t_0_err    = parameters[5]
	A_err      = parameters[6]
	D_err      = parameters[7]

	
	ans    = ((-1.0*np.exp(-1.0*(t-t_0)/tau))*(A_err))**2
	ans   += ((-1.0*A*np.exp(-1.0*(t-t_0)/tau)*(t-t_0)/(tau**2))*(tau_err))**2
	ans   += ((-1.0*A*np.exp(-1.0*(t-t_0)/tau)/tau)*(t_0_err))**2
	ans   += ((1.0*A*np.exp(-1.0*(t-t_0)/tau)/tau)*(t_err))**2
	ans    = np.sqrt(ans)

	return ans

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
	
	plt.rc('xtick',labelsize=10)
	plt.rc('ytick',labelsize=10)
	

	#Experiment
	fig = plt.figure(fig_num)
	fullPlot = fig.add_subplot(2,2,1)
	plt.plot(xdata, ydata, 'ko', markersize=1, color='#3F7F4C', label="Simulation Data")
	plt.fill_between(xdata, ydata-yerror, ydata+yerror, alpha=.8, edgecolor='#3F7F4C', facecolor='#7EFF99', linewidth=0)

	#Best Theory
	theory_label = "Theory:"

	
	theory_label += " " + fit_eq

	for i in range(len(params)):
		theory_label += "\n"
		theory_label += params_names[i] + ": "
		theory_label += ("{:.2E}".format(Decimal(params[i])))
		if str(params_err[i]) != "":
			theory_label += " $\\pm$ " + (("{:.2E}".format(Decimal(params_err[i]))))

	plt.plot(xtheory, ytheory, 'k-', color='#CC4F1B', label=theory_label)
	
	try:
		if (graph_labels['custom_placement']):
			plt.legend(loc=graph_labels['custom_placement'], prop={'size': 8})
	except:
		plt.legend(loc=1, prop={'size': 8.5})
	plt.xlabel(x_label, fontsize=12)
	plt.ylabel(y_label, fontsize=12)
	plt.title(title, fontsize=14)

	return fig_num
	


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
	plt.rc('xtick',labelsize=12)
	plt.rc('ytick',labelsize=12)

	print("xdata len: " + str(len(xdata)) + " " + "ydata len: " + str(len(ydata)))
	# print("residuals(data): " + str(ydata))
	# print("residuals(theory): " + str(ytheory))
	yres    = [(ydata[i] - ytheory[i])/yerror[i] for i in range(data_points)]
	xres    = xdata
	fig     = plt.figure(fig_num)
	resPlot = fig.add_subplot(2,2,3)
	
	plt.plot(xres, yres, 'k-', color='#3F7F4C')
	plt.xlabel(x_label, fontsize=12)
	plt.ylabel(y_label, fontsize=12)
	plt.title(title + " (Residuals)", fontsize=14)

	# #INSET
	# insetPlot = fig.add_subplot(2,2,3)
	


	# plt.xlabel(x_label, fontsize=12)
	# plt.ylabel(y_label, fontsize=12)
	# plt.title(title + " (Inset)", fontsize=14)

	# try:
	# 	if(graph_labels['reduce_temp']): #note: will be false if reduce_temp = 0
	# 		T_c   = graph_labels['reduce_temp']
	# 		xdata = [np.abs((T - T_c)/T_c) for T in xdata]
	# 		xtheory = [np.abs((T - T_c)/T_c) for T in xtheory]
	# 		plt.xlabel("Reduced Temperature $|t| = \\frac{|T - T_c|}{T_c}$ (Unitless)", fontsize=12)
	# except:
	# 	pass

	# try:
	# 	if(graph_labels['log_scale_x'] and graph_labels['log_scale_y']):
	# 			plt.plot(xdata, ydata, 'ko', markersize=1, color='#3F7F4C', label="Simulation Data")
	# 			plt.plot(xtheory, ytheory, 'k-', color='#CC4F1B')
	# 			plt.fill_between(xdata, ydata-yerror, ydata+yerror, alpha=.8, edgecolor='#3F7F4C', facecolor='#7EFF99', linewidth=0)
	# 			plt.xscale("log")
	# 			plt.yscale("log")
	# 			print("LOG INSET " + title)

	# except:
	# 	plt.plot(xdata, ydata, 'ko', markersize=1, color='#3F7F4C', label="Simulation Data")
	# 	plt.plot(xtheory, ytheory, 'k-', color='#CC4F1B')
	# 	plt.fill_between(xdata, ydata-yerror, ydata+yerror, alpha=.8, edgecolor='#3F7F4C', facecolor='#7EFF99', linewidth=0)
	
	
	
	title = title.replace("$", "")
	title = title.replace("/", ":")
	

	if (savefigs):
		fig = plt.gcf()
		fig.set_size_inches((15, 9.5), forward=False)
		fig.savefig("figs/" + title + ".png", dpi=500)

	return fig_num  + 1

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
#-----------------------------------------------------------------------------------------------------------------------------------------


#Fitting
#-----------------------------------------------------------------------------------------------------------------------------------------
fig_num                   = 1

#-----------------------------------------------------------------------------------------------------------------------------------------
T_2_set               = []
T_2_err_set           = []
concentration_set     = []



for iteration in range(num_files):
	#PI PULSE FITTING (t < t_pulse):
	concentration   = concentrations[iteration] #concentration
	print("concentration = " + str(concentration) + "\n")
	print("\tCPMG Fit : BEGIN\n")

	Estimate_tau    = 10
	Estimate_t_0    = .018
	Estimate_A      = 1
	Estimate_D      = 1

	parameter_estimates       = [Estimate_tau, Estimate_t_0, Estimate_A, Estimate_D]
	var                       = [[-np.inf, np.inf] for x in range(len(parameter_estimates))]
	parameter_bound_ranges    = ([0, -np.inf, 0, -np.inf], [np.inf, np.inf, np.inf, np.inf])

	print(parameter_bound_ranges)
	# #FIT X DATA (t < t_sweep)
	# #--------------------------------------------------------------------------------------------------------------------------------------
	xdata_data                = times[iteration]
	ydata_data                = sigma_Y[iteration]
	yerror_data               = sigma_Y_err[iteration]
	data_points               = len(xdata_data)
	function                  = CPMG_func
	


	#TAKE RANGE OF TOTAL DATA TO FIT
	xdata_fit                 = xdata_data
	ydata_fit                 = ydata_data
	yerror_fit                = yerror_data

	# plt.plot(xdata_fit, function(xdata_fit, *parameter_estimates), 'g')
	# plt.errorbar(xdata_data, ydata_data,yerror_data, fmt='b')
	# plt.errorbar(xdata_fit, ydata_fit,yerror_fit, fmt='r')
	# plt.show()
	
	popt, perr            = optimal_fit(function, xdata_fit, ydata_fit, yerror_fit, parameter_estimates, parameter_bound_ranges)
	Optimal_params        = [x for x in popt]
	Optimal_params_err    = [x for x in perr]
	Optimal_params_names  = ["$\\tau$", "$t_0$", "$A$", "$D$"]


	chi_square_min        = chi_square(function, xdata_fit, ydata_fit, yerror_fit, *Optimal_params)
	Optimal_params.append(chi_square_min)
	Optimal_params_err.append("")
	Optimal_params_names.append("Min Chi Square")


	xfit = xdata_fit
	yfit = function(xfit, *Optimal_params)

	title    = "$\\sigma_y(t)$ vs Time for concentration = " + str(concentration) + "%"
	y_label  = "$\\sigma_y(t)$ (V)"
	x_label  = "Time $t$ (s)"
	fit_eq   = "$-Ae^{-(t - t_0)/\\tau} + D$"

	fig_num = plot_fit(xdata_data, ydata_data, yerror_data, xfit, yfit, Optimal_params, Optimal_params_err, Optimal_params_names, fig_num, title=title, x_label=x_label, y_label=y_label, fit_eq=fit_eq, custom_placement=4)
	fig_num = plot_residuals(xdata_fit, ydata_fit, yerror_fit, xfit, yfit, fig_num, title=title, x_label=x_label, y_label=y_label)



	#SAVE RUN DATA
	T_2_set.append(Optimal_params[0])        
	T_2_err_set.append(Optimal_params_err[0])
	concentration_set.append(concentration)
	plt.show()
	print("\tCPMG Fit: END\n")
plt.figure(10)
plt.xlabel("Concentration (%)", fontsize=12)
plt.ylabel("$\\tau$ (s)", fontsize=12)
plt.title( "$\\tau$ vs Concentration", fontsize=14)
plt.plot(concentration_set, T_2_set)
if (savefigs):
		fig = plt.gcf()
		fig.set_size_inches((15.5, 8.5), forward=False)
		fig.savefig("figs/" + "tau vs Concentration" + ".png", dpi=500)
plt.show()
#-----------------------------------------------------------------------------------------------------------------------------------------