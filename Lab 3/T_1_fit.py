import numpy as np
import scipy.optimize as opt
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import math
import csv
import scipy.special as special
import scipy.signal as signal
from random import randint


# For displaying values
from decimal import Decimal

#DATA PARSING
#-----------------------------------------------------------------------------------------------------------------------------------------
#T,E,M,C,X  - Temperature, Energy, Magnetism, Specific Heat, Susceptibility
current_date  = '20190424'
savedata = True
savefigs = True

# concentration = .50
# name          = 'Data/20190416/t1/0416_s5_t1_'
# time_pulses = [.029, .058, .087, .116, .145, .174]
# start_point = .0007

# concentration = .60
# name          = 'Data/20190416/t1/0416_s6_t1_'
# time_pulses = [.022, .045, .067, .090, .112, .134]
# start_point = .0007

# concentration = .70
# name          = 'Data/20190416/t1/0416_s7_t1_'
# time_pulses = [.017, .034, .051, .068, .085, .102]
# start_point = .00069


# concentration = .80
# name          = 'Data/20190416/t1/0416_s8_t1_'
# time_pulses = [.014, .028, .042, .056, .070, .084]
# start_point = .0006

concentration = .90
name          = 'Data/20190416/t1/0416_s9_t1_'
time_pulses = [.008, .016, .024, .0335, .040, .048]
start_point = .00045
datafiles   = [name + str(int(1000*x)) + 'ms.csv' for x in time_pulses]
datafiles[3]= name+'33_5'+'ms.csv'

# datafiles   = [name + str(int(1000*x)) + 'ms.csv' for x in time_pulses]

print(datafiles)
num_files   = len(datafiles)

times     = []
sigma_X   = []
sigma_Y   = []
sigma     = []
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
	sigma_i     = [np.sqrt(float(row[1])**(2.0) + float(row[2])**(2.0)) for row in reader]

	times_i     = np.array(times_i).astype(np.float)
	sigma_X_i   = np.array(sigma_X_i).astype(np.float)
	sigma_Y_i   = np.array(sigma_Y_i).astype(np.float)
	sigma_i     = np.array(sigma_i).astype(np.float)

	times.append(times_i)
	sigma_X.append(sigma_X_i)
	sigma_Y.append(sigma_Y_i)
	sigma.append(sigma_i)

times     = np.array(times)
sigma_X   = np.array(sigma_X)
sigma_Y   = np.array(sigma_Y)
sigma     = np.array(sigma)

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
	# popt, pcov = opt.curve_fit(function, xdata, ydata,  p0=parameter_estimates,  bounds=parameter_ranges, maxfev=10**7)
	popt, pcov = opt.curve_fit(function, xdata, ydata, p0=parameter_estimates, maxfev=10**8)
	perr = np.sqrt(np.diag(pcov))
	return (popt, perr)

#FIT EQUATIONS
#-----------------------------------------------------------------------------------------------------------------------------------------
#PI PULSE (T_2^*)

def offset_func(x, *parameters):
	return parameters[0]

def exp_func(x, *parameters):
	#tau, delta, t_0, A, phi, D
	t      = x
	tau    = parameters[0] 
	A      = parameters[1]
	D      = parameters[2]

	ans    = -1.0*A*np.e**(-1.0*(t)/tau) + D
	return ans

def exp_err_func(x, xerr, *parameters): #the same as pi_pulse_err_func cause t = t_0
	t      = x
	tau    = parameters[0] 
	A      = parameters[1]
	D      = parameters[2]

	t_err      = xerr
	tau_err    = parameters[4]
	A_err      = parameters[5]
	D_err      = parameters[6]

	
	ans    = ((1.0*np.exp(-1.0*(t)/tau))*(A_err))**2
	ans   += ((1.0*A*np.exp(-1.0*(t)/tau)*(t)/(tau**2))*(tau_err))**2
	ans   += ((-1.0*A*np.exp(-1.0*(t)/tau)/tau)*(t_err))**2
	ans    = np.sqrt(ans)

	return ans

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


#Account for Offset
#-----------------------------------------------------------------------------------------------------------------------------------------
for iteration in range(num_files):
	time_pulse = time_pulses[iteration]
	
	Estimate_offset = 0.0
	parameter_estimates = [Estimate_offset]
	var                       = [[-np.inf, np.inf] for x in range(len(parameter_estimates))]
	parameter_bound_ranges    = ([parameter_estimates[i] + var[i][0] for i in range(len(var))], [parameter_estimates[i] + var[i][1] for i in range(len(var))])

	#SIGMAX OFFSET
	xdata_data  = times[iteration]
	ydata_data  = sigma_X[iteration]
	data_points = len(xdata_data)
	yerror_data = [.005 for x in range(data_points)] #DAQ Error
	function    = offset_func

	#TAKE RANGE OF TOTAL DATA TO FIT
	first_time                = xdata_data[0]
	end_index                 = find_nearest(xdata_data, first_time+.002)
	xdata_fit                 = xdata_data[:end_index]
	ydata_fit                 = ydata_data[:end_index]
	yerror_fit                = yerror_data[:end_index]

	popt, perr            = optimal_fit(function, xdata_fit, ydata_fit, yerror_fit, parameter_estimates, parameter_bound_ranges)
	Optimal_params        = [x for x in popt]
	Optimal_params_err    = [x for x in perr]

	sigma_X[iteration] = sigma_X[iteration] - Optimal_params[0]

	#SIGMAY OFFSET
	xdata_data  = times[iteration]
	ydata_data  = sigma_Y[iteration]
	data_points = len(xdata_data)
	yerror_data = [.005 for x in range(data_points)] #DAQ Error
	function    = offset_func

	#TAKE RANGE OF TOTAL DATA TO FIT
	start_index               = find_nearest(xdata_data, .0212)
	xdata_fit                 = xdata_data[start_index:]
	ydata_fit                 = ydata_data[start_index:]
	yerror_fit                = yerror_data[start_index:]

	popt, perr            = optimal_fit(function, xdata_fit, ydata_fit, yerror_fit, parameter_estimates, parameter_bound_ranges)
	Optimal_params        = [x for x in popt]
	Optimal_params_err    = [x for x in perr]

	sigma_Y[iteration] = sigma_Y[iteration] - Optimal_params[0]


	sigma[iteration]   = np.sqrt(sigma_Y[iteration]**(2.0) + sigma_X[iteration]**(2.0))


#-----------------------------------------------------------------------------------------------------------------------------------------

# new_times = []
# new_sigma = []
# #Take just the peaks!
# #-----------------------------------------------------------------------------------------------------------------------------------------
# for iteration in range(num_files):
# 	# time_pulse = time_pulses[iteration]
# 	xdata_data = times[iteration]
# 	ydata_data = sigma[iteration]
# 	data_points               = len(xdata_data)
# 	yerror_data               = [.005 for x in range(data_points)]

# 	plt.errorbar(xdata_data, ydata_data, yerror_data, fmt='g')
# 	plt.show()

# 	peaks, _ = signal.find_peaks(ydata_data, height=0, prominence=.02, width=2)
# 	new_times.append(xdata_data[peaks])
# 	new_sigma.append(ydata_data[peaks])
# times = np.array(new_times)
# sigma = np.array(new_sigma)


#-----------------------------------------------------------------------------------------------------------------------------------------


# #Store T_2* data
# #-----------------------------------------------------------------------------------------------------------------------------------------
# times_spread  = np.array([])
# peaks_spread  = np.array([])

# for iteration in range(num_files):
# 	#PI PULSE FITTING (t < t_pulse):
# 	time_pulse                = time_pulses[iteration]

# 	xdata_data                = times[iteration]
# 	ydata_data                = sigma[iteration]
# 	data_points               = len(xdata_data)
# 	yerror_data               = [.005 for x in range(data_points)]
# 	function                  = exp_func

# 	#TAKE RANGE OF TOTAL DATA TO FIT
# 	start_index               = find_nearest(xdata_data, .0008)
# 	end_index                 = find_nearest(xdata_data, time_pulse-.0001)

# 	print(time_pulse - .0001)
# 	xdata_fit                 = xdata_data[start_index:end_index]
# 	ydata_fit                 = ydata_data[start_index:end_index]
# 	yerror_fit                = yerror_data[start_index:end_index]


# 	times_spread = np.concatenate((times_spread, xdata_fit), axis=0)
# 	peaks_spread = np.concatenate((peaks_spread, ydata_fit), axis=0)
# 	plt.errorbar(xdata_data, ydata_data, yerror_data, fmt='b')
# 	plt.errorbar(xdata_fit, ydata_fit, yerror_fit, fmt='r')
# 	plt.show()
# #SORT by times
# sorted_time_indices = np.argsort(times_spread)
# times_spread = times_spread[sorted_time_indices]
# peaks_spread = peaks_spread[sorted_time_indices]

# plt.show()
# #-----------------------------------------------------------------------------------------------------------------------------------------

# fig_num = 1
# #Fit T_2* Spreading
# #-----------------------------------------------------------------------------------------------------------------------------------------
# Estimate_tau    = 10.0
# Estimate_A      = 1.0
# Estimate_D      = 1.0

# parameter_estimates       = [Estimate_tau,  Estimate_A, Estimate_D]
# var                       = [[0, np.inf], [0,np.inf], [-np.inf, np.inf]]
# parameter_bound_ranges    = ([parameter_estimates[i] + var[i][0] for i in range(len(var))], [parameter_estimates[i] + var[i][1] for i in range(len(var))])

# xdata_data                = times_spread
# ydata_data                = peaks_spread
# data_points               = len(xdata_data)
# yerror_data               = [.005 for x in range(data_points)] #DAQ Error
# function                  = exp_func



# #TAKE RANGE OF TOTAL DATA TO FIT
# start_index               = 0

# xdata_fit                 = xdata_data[start_index:]
# ydata_fit                 = ydata_data[start_index:]
# yerror_fit                = yerror_data[start_index:]

# popt, perr            = optimal_fit(function, xdata_fit, ydata_fit, yerror_fit, parameter_estimates, parameter_bound_ranges)
# Optimal_params        = [x for x in popt]
# Optimal_params_err    = [x for x in perr]
# Optimal_params_names  = ["$\\tau$", "$A$", "$D$"]


# chi_square_min        = chi_square(function, xdata_fit, ydata_fit, yerror_fit, *Optimal_params)
# Optimal_params.append(chi_square_min)
# Optimal_params_err.append("")
# Optimal_params_names.append("Min Chi Square")


# xfit = xdata_fit
# yfit = function(xfit, *Optimal_params)

# title    = "$\\sigma$ vs Time after $\\pi$ pulse"
# y_label  = "$\\sigma(t) = \\sqrt{\\sigma_x(t)^2 + \\sigma_y(t)^2}$ (V)"
# x_label  = "Time $t$ (s)"
# # fit_eq   = "$\\hat{E}(T) = A*E(T) + B$"
# fit_eq   = "$-Ae^{-t/\\tau} + D$"

# fig_num = plot_fit(xdata_data, ydata_data, yerror_data, xfit, yfit, Optimal_params, Optimal_params_err, Optimal_params_names, fig_num, title=title, x_label=x_label, y_label=y_label, fit_eq=fit_eq, custom_placement=1)
# fig_num = plot_residuals(xdata_fit, ydata_fit, yerror_fit, xfit, yfit, fig_num, title=title, x_label=x_label, y_label=y_label)

# plt.show()
#-----------------------------------------------------------------------------------------------------------------------------------------

# Fit Echo Peaks
#-----------------------------------------------------------------------------------------------------------------------------------------
times_T_2 = []
peaks_T_2 = []
#ADD FIRST PEAK
# times_T_2.append(0.0)
# peaks_T_2.append(exp_func(0.0, *Optimal_params))
colors    = []
for i in range(20):
    colors.append('#%06X' % randint(0, 0xFFFFFF))

plt.figure(1)
plt.axvline(x=0.0)
for iteration in range(num_files):
	#PI PULSE FITTING (t < t_pulse):
	time_pulse                = time_pulses[iteration]

	#--------------------------------------------------------------------------------------------------------------------------------------
	xdata_data                = times[iteration]
	ydata_data                = sigma_Y[iteration]
	data_points               = len(xdata_data)
	yerror_data               = [.05 for x in range(data_points)] #DAQ Error
	function                  = pi_2_pulse_func
	

	#TAKE RANGE OF TOTAL DATA TO FIT
	start_index               = find_nearest(xdata_data, time_pulse + start_point)
	end_index                 = find_nearest(xdata_data, time_pulse + .005)

	xdata_fit                 = xdata_data[start_index:]
	ydata_fit                 = ydata_data[start_index:]
	yerror_fit                = yerror_data[start_index:]


	# plt.errorbar(xdata_data, ydata_data, yerror_data, fmt=colors[iteration*2])
	# plt.errorbar(xdata_fit, ydata_fit, yerror_fit, fmt=colors[iteration*2 + 1])
	plt.axvline(x=time_pulse)
	plt.ylim(-2.5, 4.5)
	plt.xlim(-.002, .06)
	# plt.xlim(-.002, time_pulse+.02)
	plt.xlabel("Times $t$ (s)", fontsize=12)
	plt.ylabel("$\\sigma_Y$ (V)", fontsize=12)
	plt.title( "$\\sigma_Y$ vs Time (s) ($T_1$ fit) Concentration = " + str(int(100*concentration)) + "%", fontsize=14)
	plt.errorbar(xdata_data, ydata_data, yerror_data, fmt='b')
	plt.errorbar(xdata_fit, ydata_fit, yerror_fit, fmt='g')

	ydata_abs_fit  = [np.abs(y) for y in ydata_fit]
	peak_ind = np.argmax(ydata_abs_fit)
	times_T_2.append(xdata_fit[peak_ind])
	peaks_T_2.append(ydata_fit[peak_ind])

times_T_2 = np.array(times_T_2)
peaks_T_2 = np.array(peaks_T_2)
peaks_err_T_2 = [.05 for x in range(len(peaks_T_2))]
plt.errorbar(times_T_2, peaks_T_2, peaks_err_T_2, fmt='r')
if (savefigs):
		fig = plt.gcf()
		fig.savefig("figs/" + "pulse_sequence_T_1_concentration_" + str(int(100*concentration)) + ".png", dpi=500)
plt.show()
#-----------------------------------------------------------------------------------------------------------------------------------------

fig_num=2
#Fit T_2 DECOHERENCE
#-----------------------------------------------------------------------------------------------------------------------------------------
Estimate_tau    = 10.0
Estimate_A      = 1.0
Estimate_D      = 1.0

parameter_estimates       = [Estimate_tau,  Estimate_A, Estimate_D]
var                       = [[0, np.inf], [0,np.inf], [-np.inf, np.inf]]
parameter_bound_ranges    = ([parameter_estimates[i] + var[i][0] for i in range(len(var))], [parameter_estimates[i] + var[i][1] for i in range(len(var))])

xdata_data                = times_T_2
ydata_data                = peaks_T_2
data_points               = len(xdata_data)
yerror_data               = [.005 for x in range(data_points)] #DAQ Error
function                  = exp_func



#TAKE RANGE OF TOTAL DATA TO FIT
start_index               = 0

xdata_fit                 = xdata_data[start_index:]
ydata_fit                 = ydata_data[start_index:]
yerror_fit                = yerror_data[start_index:]

popt, perr            = optimal_fit(function, xdata_fit, ydata_fit, yerror_fit, parameter_estimates, parameter_bound_ranges)
Optimal_params        = [x for x in popt]
Optimal_params_err    = [x for x in perr]
Optimal_params_names  = ["$T_1$", "$A$", "$D$"]


chi_square_min        = chi_square(function, xdata_fit, ydata_fit, yerror_fit, *Optimal_params)
Optimal_params.append(chi_square_min)
Optimal_params_err.append("")
Optimal_params_names.append("Min Chi Square")


xfit = xdata_fit
yfit = function(xfit, *Optimal_params)

title    = "$\\sigma_Y$ Extrema vs Time Concentration = " + str(int(100*concentration)) + "%"
y_label  = "$\\sigma_Y(t)$ Extrema (V)"
x_label  = "Time $t$ (s)"
# fit_eq   = "$\\hat{E}(T) = A*E(T) + B$"
fit_eq   = "$-Ae^{-t/T_1} + D$"

fig_num = plot_fit(xdata_data, ydata_data, yerror_data, xfit, yfit, Optimal_params, Optimal_params_err, Optimal_params_names, fig_num, title=title, x_label=x_label, y_label=y_label, fit_eq=fit_eq, custom_placement=1)
fig_num = plot_residuals(xdata_fit, ydata_fit, yerror_fit, xfit, yfit, fig_num, title=title, x_label=x_label, y_label=y_label)

plt.show()
#-----------------------------------------------------------------------------------------------------------------------------------------

#SAVE DATA
#-----------------------------------------------------------------------------------------------------------------------------------------
row = [concentration, Optimal_params[0], Optimal_params_err[0]]
if savedata:
	with open('Analysis/' + current_date + '_T_1_fit.csv', 'a') as csvFile:
		writer = csv.writer(csvFile)
		writer.writerow(row)
	csvFile.close()
#-----------------------------------------------------------------------------------------------------------------------------------------
