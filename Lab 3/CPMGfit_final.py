import numpy as np
import scipy.optimize as opt
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import math
import csv
import scipy.special as special
import scipy.signal as signal


# For displaying values
from decimal import Decimal

#DATA PARSING
#-----------------------------------------------------------------------------------------------------------------------------------------
#T,E,M,C,X  - Temperature, Energy, Magnetism, Specific Heat, Susceptibility
datafiles      = ['Data/20190416/CPMG/0416_s5_CPMG_25ms.csv', 'Data/20190416/CPMG/0416_s6_CPMG_25ms.csv', 'Data/20190416/CPMG/0416_s7_CPMG_15ms.csv', 'Data/20190416/CPMG/0416_s8_CPMG_15ms.csv', 'Data/20190416/CPMG/0416_s9_CPMG_15ms.csv']
concentrations = [.50, .60, .70, .80, .90]
time_pulses    = [.025,.025,.015,.015,.015]
concentrations = [x*100 for x in concentrations]
num_files   = len(datafiles)
savedata      = True #save concentration vs. decay constant T_2 data
savefigs      = True
current_date  = '20190502'

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
	# popt, pcov = opt.curve_fit(function, xdata, ydata, p0=parameter_estimates, bounds=parameter_ranges, sigma=yerror, maxfev=10**7)
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

def CPMG_func(x, *parameters):
	#tau, delta, t_0, A, phi, D
	t      = x
	tau    = parameters[0] 
	t_0    = parameters[1]
	A      = parameters[2]
	D      = parameters[3]

	ans    = 1.0*A*np.e**(-1.0*(t-t_0)/tau) + D
	return ans
	
def CPMG_err_func(x, xerr, *parameters): #the same as pi_pulse_err_func cause t = t_0
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
	end_index               = find_nearest(xdata_data, -.0002)
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
	end_index               = find_nearest(xdata_data, -.0002)
	xdata_fit                 = xdata_data[:end_index]
	ydata_fit                 = ydata_data[:end_index]
	yerror_fit                = yerror_data[:end_index]

	popt, perr            = optimal_fit(function, xdata_fit, ydata_fit, yerror_fit, parameter_estimates, parameter_bound_ranges)
	Optimal_params        = [x for x in popt]
	Optimal_params_err    = [x for x in perr]

	sigma_Y[iteration] = sigma_Y[iteration] - Optimal_params[0]


	sigma[iteration]   = np.sqrt(sigma_Y[iteration]**(2.0) + sigma_X[iteration]**(2.0))


#-----------------------------------------------------------------------------------------------------------------------------------------

new_times = []
new_sigma = []
times_spread  = [None for x in range(num_files)]
peaks_spread  = [None for x in range(num_files)]
#Take just the peaks!
#Store T_2* data
#-----------------------------------------------------------------------------------------------------------------------------------------
fig_num = 1
for iteration in range(num_files):
	concentration = concentrations[iteration]
	plt.figure(fig_num)
	fig_num += 1
	#Take just the peaks!
	time_pulse = time_pulses[iteration]
	xdata_data = times[iteration]
	ydata_data = sigma[iteration]
	data_points               = len(xdata_data)
	yerror_data               = np.array([.005 for x in range(data_points)])
	first_pi_pulse_index = find_nearest(xdata_data, time_pulse)
	if (time_pulse > .02):
		peaks, _ = signal.find_peaks(ydata_data, height=0, prominence=.02, width=10)
		peaks = peaks[1:len(peaks)-1]
	else:
		peaks, _ = signal.find_peaks(ydata_data, height=0, prominence=.01, width=15)
		if(concentrations[iteration] == 80 or concentrations[iteration] == 70):
			peaks = peaks[0:len(peaks)-1]	
		else:
			peaks = peaks[1:len(peaks)-1]
	# peaks, _ = signal.find_peaks(ydata_data[first_pi_pulse_index:], height=0, prominence=.02, width=2)
	# peaks2, _= signal.find_peaks(ydata_data[0:first_pi_pulse_index], height=0, prominence=.001, width=.1)
	# peaks = np.concatenate((peaks2, peaks), axis=0)
	
	new_times.append(xdata_data[peaks])
	new_sigma.append(ydata_data[peaks])


	#Store T_2* data
	#PI PULSE FITTING (t < t_pulse):
	time_pulse                = time_pulses[iteration]/2.0

	xdata_data                = times[iteration]
	ydata_data                = sigma[iteration]
	data_points               = len(xdata_data)
	yerror_data               = np.array([.005 for x in range(data_points)])
	function                  = exp_func

	#TAKE RANGE OF TOTAL DATA TO FIT
	if(time_pulse > .02):
		start_index               = find_nearest(xdata_data, .0015)
		end_index                 = find_nearest(xdata_data, time_pulse-.008)
	else:
		start_index               = find_nearest(xdata_data, .0015)
		end_index                 = find_nearest(xdata_data, time_pulse-.0001)

	print(time_pulse - .0001)
	xdata_fit                 = xdata_data[start_index:end_index]
	ydata_fit                 = ydata_data[start_index:end_index]
	yerror_fit                = yerror_data[start_index:end_index]


	times_spread[iteration] = xdata_fit
	peaks_spread[iteration] = ydata_fit
	plt.axvline(x=0.0, color='b')
	plt.axvline(x=time_pulse, color='b')
	
	# #Straight lines at pulses
	i = 2.0
	while time_pulse*(i + 1.0) < xdata_data[len(xdata_data)-1]:
		plt.axvline(x=time_pulse*(i + 1.0), color='b')
		i += 2.0

	if(concentrations[iteration] == 90):
		plt.xlim([0.0, .10])	
	else:
		plt.xlim([0.0, time_pulse*20])
	plt.ylim(-0.5, 6)
	plt.xlabel("Time $t$ (s)", fontsize=12)
	plt.ylabel("$\\sigma(t) = \\sqrt{\\sigma_X(t)^2 + \\sigma_Y(t)^2}$ (s)", fontsize=12)
	plt.title( "$\\sigma(t)$ vs Time (s) (CPMG fit) for Concentration = " + str(concentration) + "%", fontsize=14)
	plt.errorbar(xdata_data, ydata_data, yerror_data, fmt='b')
	plt.errorbar(xdata_fit, ydata_fit, yerror_fit, fmt='r')
	plt.errorbar(xdata_data[peaks], ydata_data[peaks], yerror_data[peaks], fmt='r')
	if (savefigs):
		fig = plt.gcf()
		fig.savefig("figs/" + "CPMG_pulse_sequence_concentration_" + str(concentration) + ".png", dpi=500)
	plt.show()

times = np.array(new_times)
sigma = np.array(new_sigma)	
	



#-----------------------------------------------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------------------------------------------
#SORT by times
#-----------------------------------------------------------------------------------------------------------------------------------------
#Switch to peaks after T_2* data collection
T_2_star = []
T_2_star_err = []
#Fit T_2* Spreading
#-----------------------------------------------------------------------------------------------------------------------------------------
for iteration in range(num_files):
	Estimate_tau    = 10.0
	Estimate_A      = 1.0
	Estimate_D      = 1.0
	concentration   = concentrations[iteration]
	parameter_estimates       = [Estimate_tau,  Estimate_A, Estimate_D]
	var                       = [[-np.inf, np.inf], [-np.inf,np.inf], [-np.inf, np.inf]]
	parameter_bound_ranges    = ([parameter_estimates[i] + var[i][0] for i in range(len(var))], [parameter_estimates[i] + var[i][1] for i in range(len(var))])

	xdata_data                = times_spread[iteration]
	ydata_data                = peaks_spread[iteration]
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
	Optimal_params_names  = ["$T_2^*$", "$A$", "$D$"]


	chi_square_min        = chi_square(function, xdata_fit, ydata_fit, yerror_fit, *Optimal_params)
	Optimal_params.append(chi_square_min)
	Optimal_params_err.append("")
	Optimal_params_names.append("Min Chi Square")


	xfit = xdata_fit
	yfit = function(xfit, *Optimal_params)

	title    = "$\\sigma$ vs Time after $\\pi/2$ pulse ($T_2^*$ Fit) for Concentration = " + str(concentration) + "%"
	y_label  = "$\\sigma(t) = \\sqrt{\\sigma_x(t)^2 + \\sigma_y(t)^2}$ (V)"
	x_label  = "Time $t$ (s)"
	# fit_eq   = "$\\hat{E}(T) = A*E(T) + B$"
	fit_eq   = "$Ae^{-t/T_2^*} + D$"

	fig_num = plot_fit(xdata_data, ydata_data, yerror_data, xfit, yfit, Optimal_params, Optimal_params_err, Optimal_params_names, fig_num, title=title, x_label=x_label, y_label=y_label, fit_eq=fit_eq, custom_placement=1)
	fig_num = plot_residuals(xdata_fit, ydata_fit, yerror_fit, xfit, yfit, fig_num, title=title, x_label=x_label, y_label=y_label)
	plt.show()
	#SAVE DATA
	T_2_star.append(Optimal_params[0])
	T_2_star_err.append(Optimal_params_err[0])

plt.show()
#-----------------------------------------------------------------------------------------------------------------------------------------

#Fitting
#-----------------------------------------------------------------------------------------------------------------------------------------
# Add Zero Points
for iteration in range(num_files):
	times[iteration] = np.append(times[iteration], 1.0)
	sigma[iteration] =np.append(sigma[iteration], 0.0)
	

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
	Estimate_A      = 1
	Estimate_D      = 1

	parameter_estimates       = [Estimate_tau, Estimate_A, Estimate_D]
	var                       = [[-np.inf, np.inf] for x in range(len(parameter_estimates))]
	parameter_bound_ranges    = ([0, -np.inf, -np.inf], [np.inf, 0, np.inf])

	print(parameter_bound_ranges)
	# #FIT X DATA (t < t_sweep)
	# #--------------------------------------------------------------------------------------------------------------------------------------
	xdata_data                = times[iteration]
	ydata_data                = sigma[iteration]
	data_points               = len(xdata_data)
	yerror_data               = [.005 for x in range(data_points)] #DAQ Error
	function                  = exp_func
	


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
	Optimal_params_names  = ["$T_2$", "$A$", "$D$"]


	chi_square_min        = chi_square(function, xdata_fit, ydata_fit, yerror_fit, *Optimal_params)
	Optimal_params.append(chi_square_min)
	Optimal_params_err.append("")
	Optimal_params_names.append("Min Chi Square")


	xfit = xdata_fit
	yfit = function(xfit, *Optimal_params)

	title    = "$\\sigma(t)$ vs Time for Concentration = " + str(concentration) + "%"
	y_label  = "$\\sigma(t) = \\sqrt{\\sigma_X(t)^2 + \\sigma_Y(t)^2}$ (V)"
	x_label  = "Time $t$ (s)"
	fit_eq   = "$-Ae^{-t/T_2} + D$"

	xdata_data = xdata_data[:len(xdata_data)-1]
	ydata_data = ydata_data[:len(ydata_data)-1]
	yerror_data = yerror_data[:len(yerror_data)-1]
	xdata_fit = xdata_fit[:len(xdata_fit)-1]
	ydata_fit = ydata_fit[:len(ydata_fit)-1]
	yerror_fit = yerror_fit[:len(yerror_fit)-1]
	xfit       = xfit[:len(xfit)-1]
	yfit       = yfit[:len(yfit)-1]

	fig_num = plot_fit(xdata_data, ydata_data, yerror_data, xfit, yfit, Optimal_params, Optimal_params_err, Optimal_params_names, fig_num, title=title, x_label=x_label, y_label=y_label, fit_eq=fit_eq, custom_placement=1)
	fig_num = plot_residuals(xdata_fit, ydata_fit, yerror_fit, xfit, yfit, fig_num, title=title, x_label=x_label, y_label=y_label)



	#SAVE RUN DATA
	T_2_set.append(Optimal_params[0])        
	T_2_err_set.append(Optimal_params_err[0])
	concentration_set.append(1.0*concentration*.01)
	plt.show()
	print("\tCPMG Fit: END\n")

#-----------------------------------------------------------------------------------------------------------------------------------------

#SAVE DATA
#-----------------------------------------------------------------------------------------------------------------------------------------
rows = [[concentration_set[i], T_2_set[i], T_2_err_set[i]] for i in range(num_files)]
if savedata:
	for row in rows:
		with open('Analysis/' + current_date + '_CPMG_fit.csv', 'a') as csvFile:
			writer = csv.writer(csvFile)
			writer.writerow(row)
		csvFile.close()

rows = [[concentration_set[i], T_2_star[i], T_2_star_err[i]] for i in range(num_files)]
if savedata:
	for row in rows:
		with open('Analysis/' + current_date + '_CPMG_fit_T_2_star.csv', 'a') as csvFile:
			writer = csv.writer(csvFile)
			writer.writerow(row)
		csvFile.close()
#-----------------------------------------------------------------------------------------------------------------------------------------