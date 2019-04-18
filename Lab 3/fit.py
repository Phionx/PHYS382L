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
datafiles   = ['Data/APR:4:2019/0404_9sec.csv']
time_pulses = [.003]
num_files   = len(datafiles)

times     = []
sigma_X   = []
sigma_Y   = []
for i in range(num_files):
	datafile   = datafiles[i]
	with open(datafile, 'r') as csvFile:
		reader = csv.reader(csvFile)
		reader = list(reader)
	csvFile.close()

	time_pulse  = time_pulses[i] # time between pi/2 and pi pulse

	times_i     = [row[0] for row in reader]
	sigma_X_i   = [row[1] for row in reader]
	sigma_Y_i   = [row[2] for row in reader]

	times_i     = np.array(times_i).astype(np.float)
	sigma_X_i   = np.array(sigma_X_i).astype(np.float)
	sigma_Y_i   = np.array(sigma_Y_i).astype(np.float)
	

	times.append(times_i)
	sigma_X.append(sigma_X_i)
	sigma_Y.append(sigma_Y_i)

times     = np.array(times)
sigma_X   = np.array(sigma_X)
sigma_Y   = np.array(sigma_Y)

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
	popt, pcov = opt.curve_fit(function, xdata, ydata, p0=parameter_estimates, sigma=yerror, maxfev=10**7)
	perr = np.sqrt(np.diag(pcov))
	return (popt, perr)

#FIT EQUATIONS
#-----------------------------------------------------------------------------------------------------------------------------------------
#PI PULSE (T_2^*)

def pi_pulse_func(x, *parameters):
	#tau, delta, t_0, A, phi, D
	t      = x
	tau    = parameters[0]
	delta  = parameters[1] 
	t_0    = parameters[2]
	A      = parameters[3]
	phi    = parameters[4]
	D      = parameters[5]
	print(str([tau, delta, t_0, A, phi, D]))

	ans    = 1.0*A*np.exp(-1.0*(t-t_0)/tau)*np.sin(delta*(t-t_0) + phi) + D
	return ans

def pi_pulse_err_func(x, xerr, *parameters):
	t      = x
	tau    = parameters[0]
	Delta  = parameters[1] 
	t_0    = parameters[2]
	A      = parameters[3]
	phi    = parameters[4]
	D      = parameters[5]

	t_err      = xerr
	tau_err    = parameters[6]
	Delta_err  = parameters[7] 
	t_0_err    = parameters[8]
	A_err      = parameters[9]
	phi_err    = parameters[10]
	D_err      = parameters[11]

	ans    = ((A*Delta*np.exp(-1.0*(t-t_0)/tau)*np.cos(Delta*(t-t_0) + phi) - (A/tau)*np.exp(-1.0*(t-t_0)/tau)*np.sin(Delta*(t-t_0) + phi))**2)*(t_err)**2
	ans    += ((A/tau**2*(t-t_0)*np.exp(-1.0*(t-t_0)/tau)*np.sin(Delta*(t-t_0) + phi))**2)*(tau_err)**2
	ans    += ((A*(t-t_0)*np.exp(-(t-t_0)/tau)*np.cos(Delta*(t-t_0) + phi))**2)*(Delta_err)**2
	ans    += ((-1.0*A*Delta*np.exp(-1.0*(t-t_0)/tau)*np.cos(Delta*(t-t_0) + phi) + (A/tau)*np.exp(-1.0*(t-t_0)/tau)*np.sin(Delta*(t-t_0) + phi))**2)*(t_0_err)**2
	ans    += ((np.exp(-1.0*(t-t_0)/tau)*np.sin(Delta*(t-t_0) + phi))**2)*(A_err)**2
	ans    += ((A*np.exp(-1.0*(t-t_0)/tau)*np.cos(Delta*(t-t_0) + phi))**2)*(phi_err)**2
	ans    += ((1)**2)*(D_err)**2

	ans    = np.sqrt(ans)
	return ans


def pi_2_pulse_func(x, *parameters):
	#tau, delta, t_0, A, phi, D
	t      = x
	tau    = parameters[0]
	delta  = parameters[1] 
	t_0    = parameters[2]
	A      = parameters[3]
	phi    = parameters[4]
	D      = parameters[5]

	if (np.isscalar(t)):
		if (t < t_0):
			ans    = 1.0*A*np.exp( 1.0*(t-t_0)/tau)*np.sin(delta*(t-t_0) + phi) + D
		elif (t >= t_0):
			ans    = 1.0*A*np.exp(-1.0*(t-t_0)/tau)*np.sin(delta*(t-t_0) + phi) + D
		return ans

	T_set = x
	answer   = [0 for T in T_set]
	data_points = len(answer)
	for i in range(data_points):
		t      = T_set[i]
		if (t < t_0):
			ans    = 1.0*A*np.exp( 1.0*(t-t_0)/tau)*np.sin(delta*(t-t_0) + phi) + D
		elif (t >= t_0):
			ans    = 1.0*A*np.exp(-1.0*(t-t_0)/tau)*np.sin(delta*(t-t_0) + phi) + D
		answer[i] = ans
	answer = np.array(answer)
	return answer

	
def pi_2_pulse_err_func(x, xerr, *parameters): #the same as pi_pulse_err_func cause t = t_0
	t      = x
	tau    = parameters[0]
	Delta  = parameters[1] 
	t_0    = parameters[2]
	A      = parameters[3]
	phi    = parameters[4]
	D      = parameters[5]

	t_err      = xerr
	tau_err    = parameters[6]
	Delta_err  = parameters[7] 
	t_0_err    = parameters[8]
	A_err      = parameters[9]
	phi_err    = parameters[10]
	D_err      = parameters[11]

	
	ans    =  ((A*Delta*np.exp(-1.0*(t-t_0)/tau)*np.cos(Delta*(t-t_0) + phi) - (A/tau)*np.exp(-1.0*(t-t_0)/tau)*np.sin(Delta*(t-t_0) + phi))**2)*(t_err)**2
	ans    += ((A/tau**2*(t-t_0)*np.exp(-1.0*(t-t_0)/tau)*np.sin(Delta*(t-t_0) + phi))**2)*(tau_err)**2
	ans    += ((A*(t-t_0)*np.exp(-(t-t_0)/tau)*np.cos(Delta*(t-t_0) + phi))**2)*(Delta_err)**2
	ans    += ((-1.0*A*Delta*np.exp(-1.0*(t-t_0)/tau)*np.cos(Delta*(t-t_0) + phi) + (A/tau)*np.exp(-1.0*(t-t_0)/tau)*np.sin(Delta*(t-t_0) + phi))**2)*(t_0_err)**2
	ans    += ((np.exp(-1.0*(t-t_0)/tau)*np.sin(Delta*(t-t_0) + phi))**2)*(A_err)**2
	ans    += ((A*np.exp(-1.0*(t-t_0)/tau)*np.cos(Delta*(t-t_0) + phi))**2)*(phi_err)**2
	ans    += ((1)**2)*(D_err)**2

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
	resPlot = fig.add_subplot(2,2,2)
	
	plt.plot(xres, yres, 'k-', color='#3F7F4C')
	plt.xlabel(x_label, fontsize=12)
	plt.ylabel(y_label, fontsize=12)
	plt.title(title + " (Residuals)", fontsize=14)

	#INSET
	insetPlot = fig.add_subplot(2,2,3)
	


	plt.xlabel(x_label, fontsize=12)
	plt.ylabel(y_label, fontsize=12)
	plt.title(title + " (Inset)", fontsize=14)

	try:
		if(graph_labels['reduce_temp']): #note: will be false if reduce_temp = 0
			T_c   = graph_labels['reduce_temp']
			xdata = [np.abs((T - T_c)/T_c) for T in xdata]
			xtheory = [np.abs((T - T_c)/T_c) for T in xtheory]
			plt.xlabel("Reduced Temperature $|t| = \\frac{|T - T_c|}{T_c}$ (Unitless)", fontsize=12)
	except:
		pass

	try:
		if(graph_labels['log_scale_x'] and graph_labels['log_scale_y']):
				plt.plot(xdata, ydata, 'ko', markersize=1, color='#3F7F4C', label="Simulation Data")
				plt.plot(xtheory, ytheory, 'k-', color='#CC4F1B')
				plt.fill_between(xdata, ydata-yerror, ydata+yerror, alpha=.8, edgecolor='#3F7F4C', facecolor='#7EFF99', linewidth=0)
				plt.xscale("log")
				plt.yscale("log")
				print("LOG INSET " + title)

	except:
		plt.plot(xdata, ydata, 'ko', markersize=1, color='#3F7F4C', label="Simulation Data")
		plt.plot(xtheory, ytheory, 'k-', color='#CC4F1B')
		plt.fill_between(xdata, ydata-yerror, ydata+yerror, alpha=.8, edgecolor='#3F7F4C', facecolor='#7EFF99', linewidth=0)
	
	
	
	title = title.replace("$", "")
	title = title.replace("/", ":")
	

	if (savefigs):
		fig = plt.gcf()
		fig.set_size_inches((15.35, 9.82), forward=False)
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
X_T_2_stars           = []
X_T_2_stars_err       = []
X_T_2_stars_echo      = []
X_T_2_stars_echo_err  = []

X_deltas              = []
X_deltas_err          = []
X_deltas_echo         = []
X_deltas_echo_err     = []

sigma_X_peaks           = []
sigma_X_peaks_err       = []
sigma_X_peaks_times     = []
sigma_X_peaks_times_err = []

sigma_Y_peaks       = []
sigma_Y_peaks_times = []



for iteration in range(num_files):
	#PI PULSE FITTING (t < t_pulse):
	time_pulse                = time_pulses[i]
	print("t_sweep = " + str(time_pulse) + "\n")
	print("\tSigma_X Fit : BEGIN\n")

	Estimate_tau    = 10
	Estimate_delta  = 1
	Estimate_t_0    = .018
	Estimate_A      = 1
	Estimate_phi    = 1
	Estimate_D      = 1

	parameter_estimates       = [Estimate_tau, Estimate_delta, Estimate_t_0, Estimate_A, Estimate_phi , Estimate_D]
	var                       = [[-np.inf, np.inf] for x in range(len(parameter_estimates))]
	parameter_bound_ranges    = ([parameter_estimates[i] + var[i][0] for i in range(len(var))], [parameter_estimates[i] + var[i][1] for i in range(len(var))])

	print(parameter_bound_ranges)
	# #FIT X DATA (t < t_sweep)
	# #--------------------------------------------------------------------------------------------------------------------------------------
	# xdata_data                = times[iteration]
	# ydata_data                = sigma_X[iteration]
	# data_points               = len(xdata_data)
	# yerror_data               = [.005 for x in range(data_points)] #DAQ Error
	# function                  = pi_pulse_func
	


	# #TAKE RANGE OF TOTAL DATA TO FIT
	# start_index               = find_nearest(xdata_data, .0008)
	# end_index                 = find_nearest(xdata_data, time_pulse-.0001)

	# print(time_pulse - .0001)
	# xdata_fit                 = xdata_data[start_index:end_index]
	# ydata_fit                 = ydata_data[start_index:end_index]
	# yerror_fit                = yerror_data[start_index:end_index]

	# # plt.plot(xdata_fit, function(xdata_fit, *parameter_estimates), 'g')
	# # plt.errorbar(xdata_data, ydata_data,yerror_data, fmt='b')
	# # plt.errorbar(xdata_fit, ydata_fit,yerror_fit, fmt='r')
	# # plt.show()
	
	# popt, perr            = optimal_fit(function, xdata_fit, ydata_fit, yerror_fit, parameter_estimates, parameter_bound_ranges)
	# Optimal_params        = [x for x in popt]
	# Optimal_params_err    = [x for x in perr]
	# Optimal_params_names  = ["$\\tau$", "$\\Delta$", "$t_0$", "$A$", "$\\phi$", "$D$"]


	# chi_square_min        = chi_square(function, xdata_fit, ydata_fit, yerror_fit, *Optimal_params)
	# Optimal_params.append(chi_square_min)
	# Optimal_params_err.append("")
	# Optimal_params_names.append("Min Chi Square")


	# xfit = xdata_fit
	# yfit = function(xfit, *Optimal_params)

	# title    = "$\\sigma_x$ vs Time for $t_{sweep} = " + str(time_pulse) +  "$ (s)"
	# y_label  = "$\\sigma_x$ (V)"
	# x_label  = "Time $t$ (s)"
	# # fit_eq   = "$\\hat{E}(T) = A*E(T) + B$"
	# fit_eq   = "$Ae^{-(t - t_0)/\\tau}\\sin(\\delta(t-t_0) + \\phi) + D$"

	# fig_num = plot_fit(xdata_data, ydata_data, yerror_data, xfit, yfit, Optimal_params, Optimal_params_err, Optimal_params_names, fig_num, title=title, x_label=x_label, y_label=y_label, fit_eq=fit_eq, custom_placement=4)
	# fig_num = plot_residuals(xdata_fit, ydata_fit, yerror_fit, xfit, yfit, fig_num, title=title, x_label=x_label, y_label=y_label)

	# #SAVE RUN DATA
	# X_T_2_stars.append(Optimal_params[0])	
	# X_T_2_stars_err.append(Optimal_params_err[0])

	# X_deltas.append(Optimal_params[1])
	# X_deltas_err.append(Optimal_params_err[1])

	# X_T_2_stars.append(Optimal_params[0])
	# X_T_2_stars_err.append(Optimal_params_err[0])	
	# X_deltas.append(Optimal_params[1])
	# X_deltas_err.append(Optimal_params_err[1])

	# t_0     = Optimal_params[2]
	# t_0_err = Optimal_params_err[2]
	# sigma_X_peaks_times.append(t_0)
	# sigma_X_peaks_times_err.append(t_0_err)
	# sigma_X_peaks.append(pi_pulse_func(t_0, *Optimal_params)) #(t_0, sigma_x(t_0))
	# sigma_X_peaks_err.append(pi_pulse_err_func(t_0, t_0_err, *(Optimal_params + Optimal_params_err)))

	
	# plt.show()


	#FIT X DATA (t > t_sweep)
	#--------------------------------------------------------------------------------------------------------------------------------------
	xdata_data                = times[iteration]
	ydata_data                = sigma_X[iteration]
	data_points               = len(xdata_data)
	yerror_data               = [.05 for x in range(data_points)] #DAQ Error
	function                  = pi_2_pulse_func
	

	#TAKE RANGE OF TOTAL DATA TO FIT
	start_index               = find_nearest(xdata_data, time_pulse + .012)
	# end_index                 = find_nearest(xdata_data, time_pulse-.0001)

	xdata_fit                 = xdata_data[start_index:]
	ydata_fit                 = ydata_data[start_index:]
	yerror_fit                = yerror_data[start_index:]



	popt, perr            = optimal_fit(function, xdata_fit, ydata_fit, yerror_fit, parameter_estimates, parameter_bound_ranges)
	Optimal_params        = [x for x in popt]
	Optimal_params_err    = [x for x in perr]
	Optimal_params_names  = ["$\\tau$", "$\\Delta$", "$t_0$", "$A$", "$\\phi$", "$D$"]


	chi_square_min        = chi_square(function, xdata_fit, ydata_fit, yerror_fit, *Optimal_params)
	Optimal_params.append(chi_square_min)
	Optimal_params_err.append("")
	Optimal_params_names.append("Min Chi Square")


	xfit = xdata_fit
	yfit = function(xfit, *Optimal_params)

	title    = "$\\sigma_x$ vs Time for $t_{sweep} = " + str(time_pulse) +  "$ (s)"
	y_label  = "$\\sigma_x$ (V)"
	x_label  = "Time $t$ (s)"
	# fit_eq   = "$\\hat{E}(T) = A*E(T) + B$"
	fit_eq   = "$(t < t_{s}) \\Rightarrow Ae^{-(t - t_0)/\\tau}\\sin(\\delta(t-t_0) + \\phi) + D$\n$(t > t_{s}) \\Rightarrow Ae^{(t - t_0)/\\tau}\\sin(\\delta(t-t_0) + \\phi) + D$"

	fig_num = plot_fit(xdata_data, ydata_data, yerror_data, xfit, yfit, Optimal_params, Optimal_params_err, Optimal_params_names, fig_num, title=title, x_label=x_label, y_label=y_label, fit_eq=fit_eq, custom_placement=4)
	fig_num = plot_residuals(xdata_fit, ydata_fit, yerror_fit, xfit, yfit, fig_num, title=title, x_label=x_label, y_label=y_label)


	#SAVE RUN DATA
	X_T_2_stars.append(Optimal_params[0])	
	X_T_2_stars_err.append(Optimal_params_err[0])

	X_deltas.append(Optimal_params[1])
	X_deltas_err.append(Optimal_params_err[1])

	X_T_2_stars_echo.append(Optimal_params[0])
	X_T_2_stars_echo_err.append(Optimal_params_err[0])
	X_deltas_echo.append(Optimal_params[1])
	X_deltas_echo_err.append(Optimal_params_err[1])	

	t_0     = Optimal_params[2]
	t_0_err = Optimal_params_err[2]
	sigma_X_peaks_times.append(t_0)
	sigma_X_peaks_times_err.append(t_0_err)
	sigma_X_peaks.append(pi_pulse_func(t_0, *Optimal_params)) #(t_0, sigma_x(t_0))
	sigma_X_peaks_err.append(pi_pulse_err_func(t_0, t_0_err, *(Optimal_params + Optimal_params_err)))

	print("\tSigma_X Fit: END\n")
#-----------------------------------------------------------------------------------------------------------------------------------------


# #Magnetization Full Fit
# #-----------------------------------------------------------------------------------------------------------------------------------------
# print("Magnetization Full Fit: BEGIN\n")
# Estimate_T_c_mag          = 2.269

# function                  =	mag_func
# xdata_data                = T_M
# ydata_data                = M_mean
# yerror_data               = M_std

# fit_fraction_left_outer   = 1.0
# fit_fraction_left_inner   = 1.0/30

# fit_fraction_right_outer  = 1.0
# fit_fraction_right_inner  = 1.0/58

# fit_center_index          = find_nearest(xdata_data, Estimate_T_c_mag)
# data_points               = len(xdata_data)
# xdata_fit_left            = xdata_data[fit_center_index - int(fit_fraction_left_outer*fit_center_index):fit_center_index - int(fit_fraction_left_inner*fit_center_index)]
# xdata_fit_right           = xdata_data[fit_center_index + int(fit_fraction_right_inner*(data_points-fit_center_index-1)):fit_center_index + int(fit_fraction_right_outer*(data_points-fit_center_index))]
# xdata_fit                 = np.concatenate((xdata_fit_left, xdata_fit_right), axis=None)
# ydata_fit_left            = ydata_data[fit_center_index - int(fit_fraction_left_outer*fit_center_index):fit_center_index - int(fit_fraction_left_inner*fit_center_index)]
# ydata_fit_right           = ydata_data[fit_center_index + int(fit_fraction_right_inner*(data_points-fit_center_index-1)):fit_center_index + int(fit_fraction_right_outer*(data_points-fit_center_index))]
# ydata_fit                 = np.concatenate((ydata_fit_left, ydata_fit_right), axis=None)

# yerror_fit_left            = yerror_data[fit_center_index - int(fit_fraction_left_outer*fit_center_index):fit_center_index - int(fit_fraction_left_inner*fit_center_index)]
# yerror_fit_right           = yerror_data[fit_center_index + int(fit_fraction_right_inner*(data_points-fit_center_index-1)):fit_center_index + int(fit_fraction_right_outer*(data_points-fit_center_index))]
# yerror_fit                 = np.concatenate((yerror_fit_left, yerror_fit_right), axis=None)


# Optimal_params        = []
# Optimal_params_err    = []
# Optimal_params_names  = []

# chi_square_min        = chi_square(function, xdata_fit, ydata_fit, yerror_fit)
# Optimal_params.append(chi_square_min)
# Optimal_params_err.append("")
# Optimal_params_names.append("Min Chi Square")


# xfit = xdata_fit
# yfit = [function(x) for x in xfit]


# title    = "Magnetization vs. Temperature (Analytic Curve)"
# y_label  = "Magnetization $M$ (Unitless)"
# x_label  = "Temperature $T$ ($J/k_B$)"
# fit_eq   = "$M(T)$"

# fig_num = plot_fit(xdata_data, ydata_data, yerror_data, xfit, yfit, Optimal_params, Optimal_params_err, Optimal_params_names, fig_num, title=title, x_label=x_label, y_label=y_label, fit_eq=fit_eq)
# fig_num = plot_residuals(xdata_fit, ydata_fit, yerror_fit, xfit, yfit, fig_num, title=title, x_label=x_label, y_label=y_label)
# print("Magnetization Full Fit: END\n")
# #-----------------------------------------------------------------------------------------------------------------------------------------


# #Magnetization: T_C, \beta
# #-----------------------------------------------------------------------------------------------------------------------------------------
# print("Magnetization Power Law Fit: BEGIN\n")
# Estimate_T_c_beta         = 2.269
# Estimate_beta             = 1.0/8
# Constant                  = 1.0
# parameter_estimates       = [Estimate_T_c_beta, Estimate_beta, Constant]
# var                       = [[-.2, 0, -10000],[.2, .3, 10000]]

# parameter_bound_ranges    = ([parameter_estimates[i] + var[0][i] for i in range(len(var[0]))], [parameter_estimates[i] + var[1][i] for i in range(len(var[1]))])
# function                  = beta_func
# xdata_data                = T_M
# ydata_data                = M_mean
# yerror_data               = M_std


# parameter_estimates[0]    = 2.3
# fit_center_index          = find_nearest(xdata_data, parameter_estimates[0])

# fit_fraction_left_outer   = 1.0/3
# fit_fraction_left_inner   = 1.0/10
# fit_fraction_right_outer  = 0
# fit_fraction_right_inner  = 0


# data_points               = len(xdata_data)
# xdata_fit_left            = xdata_data[fit_center_index - int(fit_fraction_left_outer*fit_center_index):fit_center_index - int(fit_fraction_left_inner*fit_center_index)]
# xdata_fit_right           = xdata_data[fit_center_index + int(fit_fraction_right_inner*(data_points-fit_center_index-1)):fit_center_index + int(fit_fraction_right_outer*(data_points-fit_center_index))]
# xdata_fit                 = np.concatenate((xdata_fit_left, xdata_fit_right), axis=None)
# ydata_fit_left            = ydata_data[fit_center_index - int(fit_fraction_left_outer*fit_center_index):fit_center_index - int(fit_fraction_left_inner*fit_center_index)]
# ydata_fit_right           = ydata_data[fit_center_index + int(fit_fraction_right_inner*(data_points-fit_center_index-1)):fit_center_index + int(fit_fraction_right_outer*(data_points-fit_center_index))]
# ydata_fit                 = np.concatenate((ydata_fit_left, ydata_fit_right), axis=None)

# yerror_fit_left           = yerror_data[fit_center_index - int(fit_fraction_left_outer*fit_center_index):fit_center_index - int(fit_fraction_left_inner*fit_center_index)]
# yerror_fit_right          = yerror_data[fit_center_index + int(fit_fraction_right_inner*(data_points-fit_center_index-1)):fit_center_index + int(fit_fraction_right_outer*(data_points-fit_center_index))]
# yerror_fit                = np.concatenate((yerror_fit_left, yerror_fit_right), axis=None)


# popt, perr            = optimal_fit(function, xdata_fit, ydata_fit, yerror_fit, parameter_estimates, parameter_bound_ranges)
# Optimal_params        = [x for x in popt]
# Optimal_params_err    = [x for x in perr]
# Optimal_params_names  = ["$T_c$", "$\\beta$", "Constant $C$"]

# chi_square_min        = chi_square(function, xdata_fit, ydata_fit, yerror_fit, *Optimal_params)
# Optimal_params.append(chi_square_min)
# Optimal_params_err.append("")
# Optimal_params_names.append("Min Chi Square")


# xfit = xdata_fit
# yfit = function(xfit, *Optimal_params)

# title    = "Magnetization vs. Temperature (Power Law Fit)"
# y_label  = "Magnetization $M$ (Unitless)"
# x_label  = "Temperature $T$ ($J/k_B$)"
# fit_eq   = "$|M| = C|t|^{\\beta},\\ t = \\frac{T - T_c}{T_c}$"

# fig_num =plot_fit(xdata_data, ydata_data, yerror_data, xfit, yfit, Optimal_params, Optimal_params_err, Optimal_params_names, fig_num, title=title, x_label=x_label, y_label=y_label, fit_eq=fit_eq)
# fig_num = plot_residuals(xdata_fit, ydata_fit, yerror_fit, xfit, yfit, fig_num, title=title, x_label=x_label, y_label=y_label, log_scale_x=True, log_scale_y=True, reduce_temp=Optimal_params[0])
# print("Magnetization Power Law Fit: END\n")
# #-----------------------------------------------------------------------------------------------------------------------------------------



# #Susceptibility: T_C, \gamma
# #-----------------------------------------------------------------------------------------------------------------------------------------
# print("Susceptibility Power Law Fit: BEGIN\n")
# Estimate_T_c_gamma        = 2.269
# Estimate_gamma            = 7.0/4
# Constant                  = 1.0
# parameter_estimates       = [Estimate_T_c_gamma, Estimate_gamma, Constant]
# var                       = [[-.5, -.5, -100000],[.5, .5, 100000]]

# parameter_bound_ranges    = ([parameter_estimates[i] + var[0][i] for i in range(len(var[0]))], [parameter_estimates[i] + var[1][i] for i in range(len(var[1]))])
# function                  = gamma_func
# xdata_data                = T_X
# ydata_data                = X_mean
# yerror_data               = X_std



# # parameter_estimates[0]    = 2.3
# fit_center_index          = find_nearest(xdata_data, parameter_estimates[0])

# fit_fraction_left_outer   = 0
# fit_fraction_left_inner   = 0
# fit_fraction_right_outer  = 1.0/3
# fit_fraction_right_inner  = 1.0/10

# data_points               = len(xdata_data)
# xdata_fit_left            = xdata_data[fit_center_index - int(fit_fraction_left_outer*fit_center_index):fit_center_index - int(fit_fraction_left_inner*fit_center_index)]
# xdata_fit_right           = xdata_data[fit_center_index + int(fit_fraction_right_inner*(data_points-fit_center_index-1)):fit_center_index + int(fit_fraction_right_outer*(data_points-fit_center_index))]
# xdata_fit                 = np.concatenate((xdata_fit_left, xdata_fit_right), axis=None)
# ydata_fit_left            = ydata_data[fit_center_index - int(fit_fraction_left_outer*fit_center_index):fit_center_index - int(fit_fraction_left_inner*fit_center_index)]
# ydata_fit_right           = ydata_data[fit_center_index + int(fit_fraction_right_inner*(data_points-fit_center_index-1)):fit_center_index + int(fit_fraction_right_outer*(data_points-fit_center_index))]
# ydata_fit                 = np.concatenate((ydata_fit_left, ydata_fit_right), axis=None)

# yerror_fit_left           = yerror_data[fit_center_index - int(fit_fraction_left_outer*fit_center_index):fit_center_index - int(fit_fraction_left_inner*fit_center_index)]
# yerror_fit_right          = yerror_data[fit_center_index + int(fit_fraction_right_inner*(data_points-fit_center_index-1)):fit_center_index + int(fit_fraction_right_outer*(data_points-fit_center_index))]
# yerror_fit                = np.concatenate((yerror_fit_left, yerror_fit_right), axis=None)


# popt, perr            = optimal_fit(function, xdata_fit, ydata_fit, yerror_fit, parameter_estimates, parameter_bound_ranges)
# Optimal_params        = [x for x in popt]
# Optimal_params_err    = [x for x in perr]
# Optimal_params_names  = ["$T_c$", "$\\gamma$", "Constant $C$"]

# chi_square_min        = chi_square(function, xdata_fit, ydata_fit, yerror_fit, *Optimal_params)
# Optimal_params.append(chi_square_min)
# Optimal_params_err.append("")
# Optimal_params_names.append("Min Chi Square")


# xfit = xdata_fit
# yfit = function(xfit, *Optimal_params)

# title    = "Susceptibility vs. Temperature (Power Law Fit)"
# y_label  = "Susceptibility $\\chi$ (1/J)"
# x_label  = "Temperature $T$ (J/$k_B$)"
# fit_eq   = "$\\chi = C|t|^{-\\gamma},\\ t = \\frac{T - T_c}{T_c}$"

# fig_num = plot_fit(xdata_data, ydata_data, yerror_data, xfit, yfit, Optimal_params, Optimal_params_err, Optimal_params_names, fig_num, title=title, x_label=x_label, y_label=y_label, fit_eq=fit_eq)
# fig_num = plot_residuals(xdata_fit, ydata_fit, yerror_fit, xfit, yfit, fig_num, title=title, x_label=x_label, y_label=y_label, log_scale_x=True, log_scale_y=True, reduce_temp=Optimal_params[0])
# print("Susceptibility Power Law Fit: END\n")
# #-----------------------------------------------------------------------------------------------------------------------------------------

# #Specific Heat: T_c, \alpha
# #-----------------------------------------------------------------------------------------------------------------------------------------
# print("Specific Heat Power Law Fit: BEGIN\n")


# title    = "Specific Heat vs. Temperature (Power Law Fit: Less than $T_c$)"
# y_label  = "Specific Heat $c_V$ ($k_B$)"
# x_label  = "Temperature $T$ ($J/k_B$)"
# fit_eq   = "$c_V = C|t|^{-\\alpha},\\ t = \\frac{T - T_c}{T_c}$"

# #LEFT SIDE
# #--------------------------------------------------------------------------------------------------
# Estimate_T_c_alpha        = 2.269
# Estimate_alpha            = 0.0
# Constant                  = 1.0
# parameter_estimates       = [Estimate_T_c_alpha, Estimate_alpha, Constant]
# var                       = [[-.5, -.5, -100000],[.5,.5, 100000]]

# parameter_bound_ranges    = ([parameter_estimates[i] + var[0][i] for i in range(len(var[0]))], [parameter_estimates[i] + var[1][i] for i in range(len(var[1]))])
# function                  = alpha_func
# xdata_data                = T_C
# ydata_data                = C_mean
# yerror_data               = C_std



# # parameter_estimates[0]    = 2.3
# fit_center_index          = find_nearest(xdata_data, parameter_estimates[0])

# fit_fraction_left_outer   = 1.0/4
# fit_fraction_left_inner   = 1.0/10
# fit_fraction_right_outer  = 0.0/3
# fit_fraction_right_inner  = 0.0/10



# data_points               = len(xdata_data)
# xdata_fit_left            = xdata_data[fit_center_index - int(fit_fraction_left_outer*fit_center_index):fit_center_index - int(fit_fraction_left_inner*fit_center_index)]
# xdata_fit_right           = xdata_data[fit_center_index + int(fit_fraction_right_inner*(data_points-fit_center_index-1)):fit_center_index + int(fit_fraction_right_outer*(data_points-fit_center_index))]
# xdata_fit                 = np.concatenate((xdata_fit_left, xdata_fit_right), axis=None)
# ydata_fit_left            = ydata_data[fit_center_index - int(fit_fraction_left_outer*fit_center_index):fit_center_index - int(fit_fraction_left_inner*fit_center_index)]
# ydata_fit_right           = ydata_data[fit_center_index + int(fit_fraction_right_inner*(data_points-fit_center_index-1)):fit_center_index + int(fit_fraction_right_outer*(data_points-fit_center_index))]
# ydata_fit                 = np.concatenate((ydata_fit_left, ydata_fit_right), axis=None)

# yerror_fit_left           = yerror_data[fit_center_index - int(fit_fraction_left_outer*fit_center_index):fit_center_index - int(fit_fraction_left_inner*fit_center_index)]
# yerror_fit_right          = yerror_data[fit_center_index + int(fit_fraction_right_inner*(data_points-fit_center_index-1)):fit_center_index + int(fit_fraction_right_outer*(data_points-fit_center_index))]
# yerror_fit                = np.concatenate((yerror_fit_left, yerror_fit_right), axis=None)


# popt, perr            = optimal_fit(function, xdata_fit, ydata_fit, yerror_fit, parameter_estimates, parameter_bound_ranges)
# Optimal_params        = [x for x in popt]
# Optimal_params_err    = [x for x in perr]
# Optimal_params_names  = ["$T_c$ (Left)", "$\\alpha$ (Left)", "Constant $C$ (Left)"]

# chi_square_min        = chi_square(function, xdata_fit, ydata_fit, yerror_fit, *Optimal_params)
# Optimal_params.append(chi_square_min)
# Optimal_params_err.append("")
# Optimal_params_names.append("Min Chi Square (Left)")


# xfit = xdata_fit
# yfit = function(xdata_fit, *Optimal_params[:3])

# fig_num = plot_fit(xdata_data, ydata_data, yerror_data, xfit, yfit, Optimal_params, Optimal_params_err, Optimal_params_names, fig_num, title=title, x_label=x_label, y_label=y_label, fit_eq=fit_eq)
# fig_num = plot_residuals(xdata_fit, ydata_fit, yerror_fit, xfit, yfit, fig_num, title=title, x_label=x_label, y_label=y_label, log_scale_x=True, log_scale_y=True, reduce_temp=Optimal_params[0])

# #RIGHT SIDE
# title    = "Specific Heat vs. Temperature (Power Law Fit: Greater than $T_c$)"
# #--------------------------------------------------------------------------------------------------
# Estimate_T_c_alpha        = 2.269
# Estimate_alpha            = 0.0
# Constant                  = 1.0
# parameter_estimates       = [Estimate_T_c_alpha, Estimate_alpha, Constant]
# var                       = [[-.5, -.5, -100000],[.5,.5, 100000]]

# parameter_bound_ranges    = ([parameter_estimates[i] + var[0][i] for i in range(len(var[0]))], [parameter_estimates[i] + var[1][i] for i in range(len(var[1]))])
# function                  = alpha_func
# xdata_data                = T_C
# ydata_data                = C_mean
# yerror_data               = C_std



# # parameter_estimates[0]    = 2.3
# fit_center_index          = find_nearest(xdata_data, parameter_estimates[0])

# fit_fraction_left_outer   = 0.0/4
# fit_fraction_left_inner   = 0.0/10
# fit_fraction_right_outer  = 1.0/4
# fit_fraction_right_inner  = 1.0/12



# data_points               = len(xdata_data)
# xdata_fit_left            = xdata_data[fit_center_index - int(fit_fraction_left_outer*fit_center_index):fit_center_index - int(fit_fraction_left_inner*fit_center_index)]
# xdata_fit_right           = xdata_data[fit_center_index + int(fit_fraction_right_inner*(data_points-fit_center_index-1)):fit_center_index + int(fit_fraction_right_outer*(data_points-fit_center_index))]
# xdata_fit2                = np.concatenate((xdata_fit_left, xdata_fit_right), axis=None)
# ydata_fit_left            = ydata_data[fit_center_index - int(fit_fraction_left_outer*fit_center_index):fit_center_index - int(fit_fraction_left_inner*fit_center_index)]
# ydata_fit_right           = ydata_data[fit_center_index + int(fit_fraction_right_inner*(data_points-fit_center_index-1)):fit_center_index + int(fit_fraction_right_outer*(data_points-fit_center_index))]
# ydata_fit2                = np.concatenate((ydata_fit_left, ydata_fit_right), axis=None)

# yerror_fit_left           = yerror_data[fit_center_index - int(fit_fraction_left_outer*fit_center_index):fit_center_index - int(fit_fraction_left_inner*fit_center_index)]
# yerror_fit_right          = yerror_data[fit_center_index + int(fit_fraction_right_inner*(data_points-fit_center_index-1)):fit_center_index + int(fit_fraction_right_outer*(data_points-fit_center_index))]
# yerror_fit2               = np.concatenate((yerror_fit_left, yerror_fit_right), axis=None)


# popt, perr            = optimal_fit(function, xdata_fit2, ydata_fit2, yerror_fit2, parameter_estimates, parameter_bound_ranges)
# for x in popt:
# 	Optimal_params.append(x)

# for x in perr:
# 	Optimal_params_err.append(x)

# Optimal_params_names.append("$T_c$ (Right)")
# Optimal_params_names.append("$\\alpha$ (Right)")
# Optimal_params_names.append("Constant $C$ (Right)")


# chi_square_min        = chi_square(function, xdata_fit2, ydata_fit2, yerror_fit2, *Optimal_params[4:])
# Optimal_params.append(chi_square_min)
# Optimal_params_err.append("")
# Optimal_params_names.append("Min Chi Square (Right)")

# xfit = xdata_fit2
# yfit = function(xdata_fit2, *Optimal_params[4:])

# fig_num = plot_fit(xdata_data, ydata_data, yerror_data, xfit, yfit, Optimal_params[4:], Optimal_params_err[4:], Optimal_params_names[4:], fig_num, title=title, x_label=x_label, y_label=y_label, fit_eq=fit_eq)
# print(Optimal_params_names[4])
# fig_num = plot_residuals(xdata_fit2, ydata_fit2, yerror_fit2, xfit, yfit, fig_num, title=title, x_label=x_label, y_label=y_label, log_scale_x=True, log_scale_y=True, reduce_temp=Optimal_params[4])
# #--------------------------------------------------------------------------------------------------


# #BOTH SIDES
# title    = "Specific Heat vs. Temperature (Power Law Fit)"
# y_label  = "Specific Heat $c_V$ ($k_B$)"
# x_label  = "Temperature $T$ (J/$k_B$)"
# fit_eq   = "$c_V = C|t|^{-\\alpha},\\ t = \\frac{T - T_c}{T_c}$"


# xfit = np.concatenate((xdata_fit, xdata_fit2), axis=None)
# yfit = np.concatenate((function(xdata_fit, *Optimal_params[:3]), function(xdata_fit2, *Optimal_params[4:])), axis=None) 

# xdata_fit  = np.concatenate((xdata_fit, xdata_fit2), axis=None)
# ydata_fit  = np.concatenate((ydata_fit, ydata_fit2), axis=None)
# yerror_fit = np.concatenate((yerror_fit, yerror_fit2), axis=None)

# fig_num = plot_fit(xdata_data, ydata_data, yerror_data, xfit, yfit, Optimal_params, Optimal_params_err, Optimal_params_names, fig_num, title=title, x_label=x_label, y_label=y_label, fit_eq=fit_eq)
# fig_num = plot_residuals(xdata_fit, ydata_fit, yerror_fit, xfit, yfit, fig_num, title=title, x_label=x_label, y_label=y_label)
# print("Specific Heat Power Law Fit: END\n")
plt.show()
#-----------------------------------------------------------------------------------------------------------------------------------------
