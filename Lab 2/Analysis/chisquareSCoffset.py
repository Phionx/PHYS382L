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
#Spin Correlations
SCfile = 'Data/APR:3:2019/CLEANED_3_30_N100_1.20T3.00_SC_v1.csv'
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
X       = []
SC_mean = []
SC_std  = []

for i in range(56, len(readerSC)):
	row = readerSC[i]
	if row[2] != "": #checks if first element is there
		T_SC.append(row[0])
		SC_mean_row = []
		SC_std_row  = []
		X_row       = []
		for i in range(d_max):
			if row[2+i] != "":
				SC_mean_row.append(row[2+i])
				SC_std_row.append(row[d_max+2+i])
				X_row.append(i+1)
		SC_mean.append(np.array(SC_mean_row).astype(np.float))
		SC_std.append(np.array(SC_std_row).astype(np.float))
		X.append(np.array(X_row).astype(np.float))
T_SC = np.array(T_SC).astype(np.float)

# for i in range(len(readerSC)):
# 	print("i value : " + str(i) +"\n")
# 	print("xdata : " + str(X[i]) + "\n")
# 	print("ydata : " + str(SC_mean[i]) + "\n")
# 	print("yerror : " + str(SC_std[i]) + "\n")
# T_SC    = np.array(T_SC).astype(np.float)
# SC_mean = np.array(SC_mean).astype(np.float)
# SC_std  = np.array(SC_std).astype(np.float)
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
#-------------------------------------------------------------------------
#SC
#R(x) = Ce^(-x/\xi)
def xi_func(x, *parameters):
	x        = x
	xi       = parameters[0]
	scale    = parameters[1]
	sadd     = parameters[2]
	sadd2    = parameters[3]

	y = scale*(np.e**(-1.0*(x + sadd2)/(xi))) + sadd #+ add
	return y 

#\xi found for T near TC
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
#For R(x) near T_C
#R(x) = |x|^{-\eta}
#Theory  Value: 1/4
def eta_func(x, *parameters):
	eta   = parameters[0]
	constant = parameters[1]

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
	
	plt.rc('xtick',labelsize=10)
	plt.rc('ytick',labelsize=10)

	#Experiment
	fig = plt.figure(fig_num)
	fullPlot = fig.add_subplot(2,2,1)
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
	plt.legend(loc=1, prop={'size': 8})
	plt.xlabel(x_label, fontsize=10)
	plt.ylabel(y_label, fontsize=10)
	plt.title(title, fontsize=12)

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
	plt.xlabel(x_label, fontsize=10)
	plt.ylabel(y_label, fontsize=10)
	plt.title(title + " (Residuals)", fontsize=12)

	#INSET
	insetPlot = fig.add_subplot(2,2,3)



	plt.xlabel(x_label, fontsize=10)
	plt.ylabel(y_label, fontsize=10)
	plt.title(title + " (Inset)", fontsize=12)

	try:
		if(graph_labels['reduce_temp']): #note: will be false if reduce_temp = 0
			T_c   = graph_labels['reduce_temp']
			xdata = [np.abs((T - T_c)/T_c) for T in xdata]
			xtheory = [np.abs((T - T_c)/T_c) for T in xtheory]
			plt.xlabel("Reduced Temperature $t = \\frac{T - T_c}{T_c}$ (Unitless)", fontsize=10)
	except:
		pass

	try:
		if (graph_labels['y_offset']):
			offset = graph_labels['y_offset']
			ydata = [y - offset for y in ydata]
			ytheory = [y - offset for y in ytheory]
			plt.ylabel("Spin Correlation $R(x)$ Offset by " +  ("{:.2E}".format(Decimal(offset))) + " (Unitless)", fontsize=10)
	except:
		pass


	try:
		if(graph_labels['log_scale_x'] and graph_labels['log_scale_y']):
			plt.plot(xdata, ydata, 'ko', markersize=1, color='#3F7F4C', label="Simulation Data")
			plt.plot(xtheory, ytheory, 'k-', color='#CC4F1B')
			plt.fill_between(xdata, ydata-yerror, ydata+yerror, alpha=.8, edgecolor='#3F7F4C', facecolor='#7EFF99', linewidth=0)
			plt.xscale("log")
			plt.yscale("log")
			print("LOG INSET BOTH " + title)

	except:
		try:
			if(graph_labels['log_scale_y']):
				plt.plot(xdata, ydata, 'ko', markersize=1, color='#3F7F4C', label="Simulation Data")
				plt.plot(xtheory, ytheory, 'k-', color='#CC4F1B')
				plt.fill_between(xdata, ydata-yerror, ydata+yerror, alpha=.8, edgecolor='#3F7F4C', facecolor='#7EFF99', linewidth=0)
				plt.yscale("log")
				print("LOG INSET Y SCALE " + title)
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
#Estimates

# Estimate_T_c_nu      = 2.2
# Estimate_nu          = 1.0

# Estimate_eta         = 1.0/4

#Spin Correlation 
#-----------------------------------------------------------------------------------------------------------------------------------------

fig_num                   = 1

# Spin Correlation Exponential Fit at Every T: \xi
#-----------------------------------------------------------------------------------------------------------------------------------------
num_temps = len(T_SC)
XI        = []
XI_err    = []


for iteration in range(num_temps):
	print("Spin Correlation Fit " + str(T_SC[iteration]) + " : BEGIN\n")

	T = T_SC[iteration]
	Estimate_xi               = 1.0
	scale_constant            = 1.0
	add_constant              = 0.0
	add_constant2             = 0.0

	# parameter_estimates       = [Estimate_xi, scale_constant, add_constant]
	# var                       = [[-10000, -10000, -100000],[10000, 10000, 100000]]

	parameter_estimates       = [Estimate_xi, scale_constant, add_constant, add_constant2]
	var                       = [[-10000, -10000, -100, -100],[10000, 10000, 100, 100]]


	parameter_bound_ranges    = ([parameter_estimates[i] + var[0][i] for i in range(len(var[0]))], [parameter_estimates[i] + var[1][i] for i in range(len(var[1]))])
	function                  = xi_func
	xdata_data                = X[iteration]
	ydata_data                = SC_mean[iteration]
	yerror_data               = SC_std[iteration]
	

	xdata_fit  = xdata_data
	ydata_fit  = ydata_data
	yerror_fit = yerror_data

	popt, perr            = optimal_fit(function, xdata_fit, ydata_fit, yerror_fit, parameter_estimates, parameter_bound_ranges)
	Optimal_params        = [x for x in popt]
	Optimal_params_err    = [x for x in perr]

	XI.append(Optimal_params[0])
	XI_err.append(Optimal_params_err[0])

	# Optimal_params_names  = ["$\\xi$", "Constant A", "Constant B"]
	Optimal_params_names  = ["$\\xi$", "Constant A", "Constant B", "Constant C"]

	chi_square_min        = chi_square(function, xdata_fit, ydata_fit, yerror_fit, *Optimal_params)
	Optimal_params.append(chi_square_min)
	Optimal_params_err.append("")
	Optimal_params_names.append("Min Chi Square")


	xfit = xdata_fit
	yfit = function(xfit, *Optimal_params)

	title    = "Spin Correlation vs. Distance"
	y_label  = "Spin Correlation $R(x)$ (Unitless)"
	x_label  = "Distance $x$ (Spacing)"
	fit_eq   = "$R(x) = A*e^{-(x + C)/\\xi} + B$"

	if (iteration == 50 or iteration == num_temps - 50):
		title    = "Spin Correlation at $T = " + str(T) + " (J/k_B)$ vs. Distance (Exponential Fit)"
		fig_num = plot_fit(xdata_data, ydata_data, yerror_data, xfit, yfit, Optimal_params, Optimal_params_err, Optimal_params_names, fig_num, title=title, x_label=x_label, y_label=y_label, fit_eq=fit_eq)
		fig_num = plot_residuals(xdata_fit, ydata_fit, yerror_fit, xfit, yfit, fig_num, title=title, x_label=x_label, y_label=y_label, log_scale_y=True, y_offset=Optimal_params[2])
	print("Spin Correlation Fit: END\n")
	# plt.show()
XI = np.array(XI)
XI_err = np.array(XI_err)
#-----------------------------------------------------------------------------------------------------------------------------------------



# Spin Correlation Length: \nu
#-----------------------------------------------------------------------------------------------------------------------------------------

print("Spin Correlation Length: BEGIN\n")
Estimate_T_c_nu           = 2.269
Estimate_nu               = 1.0
scale_constant            = 1.0

parameter_estimates       = [Estimate_T_c_nu, Estimate_nu, scale_constant]
# var                       = [[-10, -10, -100],[10, 10, 100]]
var                       = [[-10, -10, -100],[10, 10, 100]]
# parameter_estimates       = [Estimate_xi, scale_constant]
# var                       = [[-10000, -10000],[10000, 10000]]


parameter_bound_ranges    = ([parameter_estimates[i] + var[0][i] for i in range(len(var[0]))], [parameter_estimates[i] + var[1][i] for i in range(len(var[1]))])
function                  = nu_func
xdata_data                = T_SC
ydata_data                = XI
yerror_data               = XI_err

# plt.plot(xdata_data, ydata_data)
# plt.show()

# parameter_estimates[0]    = 2.3
fit_center_index          = find_nearest(xdata_data, parameter_estimates[0])

fit_fraction_left_outer   = 0.0
fit_fraction_left_inner   = 0.0
fit_fraction_right_outer  = 1.0/4
fit_fraction_right_inner  = 1.0/32

data_points               = len(xdata_data)
xdata_fit_left            = xdata_data[fit_center_index - int(fit_fraction_left_outer*fit_center_index):fit_center_index - int(fit_fraction_left_inner*fit_center_index)]
xdata_fit_right           = xdata_data[fit_center_index + int(fit_fraction_right_inner*(data_points-fit_center_index-1)):fit_center_index + int(fit_fraction_right_outer*(data_points-fit_center_index))]
xdata_fit                 = np.concatenate((xdata_fit_left, xdata_fit_right), axis=None)
ydata_fit_left            = ydata_data[fit_center_index - int(fit_fraction_left_outer*fit_center_index):fit_center_index - int(fit_fraction_left_inner*fit_center_index)]
ydata_fit_right           = ydata_data[fit_center_index + int(fit_fraction_right_inner*(data_points-fit_center_index-1)):fit_center_index + int(fit_fraction_right_outer*(data_points-fit_center_index))]
ydata_fit                 = np.concatenate((ydata_fit_left, ydata_fit_right), axis=None)

yerror_fit_left           = yerror_data[fit_center_index - int(fit_fraction_left_outer*fit_center_index):fit_center_index - int(fit_fraction_left_inner*fit_center_index)]
yerror_fit_right          = yerror_data[fit_center_index + int(fit_fraction_right_inner*(data_points-fit_center_index-1)):fit_center_index + int(fit_fraction_right_outer*(data_points-fit_center_index))]
yerror_fit                = np.concatenate((yerror_fit_left, yerror_fit_right), axis=None)


popt, perr            = optimal_fit(function, xdata_fit, ydata_fit, yerror_fit, parameter_estimates, parameter_bound_ranges)
Optimal_params        = [x for x in popt]
Optimal_params_err    = [x for x in perr]
Optimal_params_names  = ["$T_c$","$\\nu$", "Constant C"]
# Optimal_params_names  = ["$\\xi$", "Constant A"]

chi_square_min        = chi_square(function, xdata_fit, ydata_fit, yerror_fit, *Optimal_params)

Optimal_params.append(chi_square_min)
Optimal_params_err.append("")
Optimal_params_names.append("Min Chi Square")


xfit = xdata_fit
yfit = function(xfit, *Optimal_params)

title    = "Spin Correlation Length vs. Temperature"
y_label  = "Spin Correlation Length $\\xi$ (Spacing)"
x_label  = "Temperature $T$ ($J/k_B$)"
fit_eq   = "$\\xi = C*|t|^{-\\nu},\\ t = \\frac{T - T_c}{T_c}$"

fig_num = plot_fit(xdata_data, ydata_data, yerror_data, xfit, yfit, Optimal_params, Optimal_params_err, Optimal_params_names, fig_num, title=title, x_label=x_label, y_label=y_label, fit_eq=fit_eq)
fig_num = plot_residuals(xdata_fit, ydata_fit, yerror_fit, xfit, yfit, fig_num, title=title, x_label=x_label, y_label=y_label,  log_scale_x=True, log_scale_y=True, reduce_temp=Optimal_params[0])
print("Spin Correlation Length Fit: END\n")

#-----------------------------------------------------------------------------------------------------------------------------------------


# Generate Data for eta fit
#-----------------------------------------------------------------------------------------------------------------------------------------
Optimal_T_c_nu     = Optimal_params[0]
Optimal_T_c_nu_err = Optimal_params_err[0] 

# XI_max     = XI[XI.argmax()]
# XI_max_err = XI_err[XI.argmax()] 


# print("XI_max : " + str(XI_max))
# print("XI_max_err : " + str(XI_max_err))


# average_center            = XI.argmax()
average_center            = find_nearest(T_SC, Optimal_T_c_nu)
average_left              = find_nearest(T_SC, Optimal_T_c_nu - Optimal_T_c_nu_err)
average_right             = find_nearest(T_SC, Optimal_T_c_nu + Optimal_T_c_nu_err)
# average_left              = average_center - 10
# average_right             = average_center + 10

# print("average : " + str(average_left) + ", " + str(average_center) + ", " + str(average_right))
SC_mean_centered = np.transpose(np.array(SC_mean[average_left:average_right+1]))

SC_mean_avg = []
SC_mean_std = []

for row in SC_mean_centered:
	SC_mean_avg.append(np.mean(row))
	SC_mean_std.append(np.std(row))

SC_mean_avg = np.array(SC_mean_avg)
SC_mean_std = np.array(SC_mean_std)

print(SC_mean_avg)
print(SC_mean_std)
#-----------------------------------------------------------------------------------------------------------------------------------------



# Find exponential fit near T_c
#-----------------------------------------------------------------------------------------------------------------------------------------
iteration                 = average_center

print("Spin Correlation Exp Fit near T_c ($T = " + str(T_SC[iteration]) + " (J/k_B)$): BEGIN\n")

T = T_SC[iteration]
print("T_c : " + str(T))
Estimate_xi               = 1.0
scale_constant            = 1.0
add_constant              = 0.0
add_constant2             = 0.0
# parameter_estimates       = [Estimate_xi, scale_constant, add_constant]
# var                       = [[-10000, -10000, -100000],[10000, 10000, 100000]]

parameter_estimates       = [Estimate_xi, scale_constant, add_constant, add_constant2]
var                       = [[-10000, -10000, -100, -100],[10000, 10000, 100, 100]]


parameter_bound_ranges    = ([parameter_estimates[i] + var[0][i] for i in range(len(var[0]))], [parameter_estimates[i] + var[1][i] for i in range(len(var[1]))])
function                  = xi_func
xdata_data                = X[iteration]
ydata_data                = SC_mean[iteration]
yerror_data               = SC_std[iteration]

xdata_fit  = xdata_data
ydata_fit  = ydata_data
yerror_fit = yerror_data

popt, perr            = optimal_fit(function, xdata_fit, ydata_fit, yerror_fit, parameter_estimates, parameter_bound_ranges)
Optimal_params        = [x for x in popt]
Optimal_params_err    = [x for x in perr]

xi_at_T_c     = Optimal_params[0]
xi_at_T_c_err = Optimal_params_err[0]

# Optimal_params_names  = ["$\\xi$", "Constant A", "Constant B"]
Optimal_params_names  = ["$\\xi$", "Constant A", "Constant B", "Constant C"]

chi_square_min        = chi_square(function, xdata_fit, ydata_fit, yerror_fit, *Optimal_params)
Optimal_params.append(chi_square_min)
Optimal_params_err.append("")
Optimal_params_names.append("Min Chi Square")


xfit = xdata_fit
yfit = function(xfit, *Optimal_params)

title    = "Spin Correlation near T_c ($T = " + str(T_SC[iteration]) + " (J/k_B)$) vs. Distance (Exponential Fit)"
y_label  = "Spin Correlation $R(x)$ (Unitless)"
x_label  = "Distance $x$ (Spacing)"
fit_eq   = "$R(x) = A*e^{-(x + C)/\\xi} + B$"

fig_num = plot_fit(xdata_data, ydata_data, yerror_data, xfit, yfit, Optimal_params, Optimal_params_err, Optimal_params_names, fig_num, title=title, x_label=x_label, y_label=y_label, fit_eq=fit_eq)
fig_num = plot_residuals(xdata_fit, ydata_fit, yerror_fit, xfit, yfit, fig_num, title=title, x_label=x_label, y_label=y_label, log_scale_y=True, y_offset=Optimal_params[2])
print("Spin Correlation Exp Fit near T_c: END\n")
# plt.show()
#-----------------------------------------------------------------------------------------------------------------------------------------	


# Spin Correlation Power Law fit near T_c: \eta
#-----------------------------------------------------------------------------------------------------------------------------------------

print("Spin Correlation Near T_c: BEGIN\n")
Estimate_eta               = 1.0
scale_constant            = 1.0

parameter_estimates       = [Estimate_eta, scale_constant]
var                       = [[-10, -10],[10, 10]]

# parameter_estimates       = [Estimate_xi, scale_constant]
# var                       = [[-10000, -10000],[10000, 10000]]


parameter_bound_ranges    = ([parameter_estimates[i] + var[0][i] for i in range(len(var[0]))], [parameter_estimates[i] + var[1][i] for i in range(len(var[1]))])
function                  = eta_func
xdata_data                = X[0]
ydata_data                = SC_mean_avg
yerror_data               = SC_mean_std

# plt.plot(xdata_data, ydata_data)
# plt.show()

# parameter_estimates[0]    = 2.3
# fit_center_index          = find_nearest(xdata_data, XI_max)
fit_center_index          = find_nearest(xdata_data, XI[average_center])

fit_fraction_left_outer   = 1.0
fit_fraction_left_inner   = 0.0
fit_fraction_right_outer  = 0.0
fit_fraction_right_inner  = 0.0

data_points               = len(xdata_data)
xdata_fit_left            = xdata_data[fit_center_index - int(fit_fraction_left_outer*fit_center_index):fit_center_index - int(fit_fraction_left_inner*fit_center_index)]
xdata_fit_right           = xdata_data[fit_center_index + int(fit_fraction_right_inner*(data_points-fit_center_index-1)):fit_center_index + int(fit_fraction_right_outer*(data_points-fit_center_index))]
xdata_fit                 = np.concatenate((xdata_fit_left, xdata_fit_right), axis=None)
ydata_fit_left            = ydata_data[fit_center_index - int(fit_fraction_left_outer*fit_center_index):fit_center_index - int(fit_fraction_left_inner*fit_center_index)]
ydata_fit_right           = ydata_data[fit_center_index + int(fit_fraction_right_inner*(data_points-fit_center_index-1)):fit_center_index + int(fit_fraction_right_outer*(data_points-fit_center_index))]
ydata_fit                 = np.concatenate((ydata_fit_left, ydata_fit_right), axis=None)

yerror_fit_left           = yerror_data[fit_center_index - int(fit_fraction_left_outer*fit_center_index):fit_center_index - int(fit_fraction_left_inner*fit_center_index)]
yerror_fit_right          = yerror_data[fit_center_index + int(fit_fraction_right_inner*(data_points-fit_center_index-1)):fit_center_index + int(fit_fraction_right_outer*(data_points-fit_center_index))]
yerror_fit                = np.concatenate((yerror_fit_left, yerror_fit_right), axis=None)

print(str(xdata_fit))
print(str(ydata_fit))

popt, perr            = optimal_fit(function, xdata_fit, ydata_fit, yerror_fit, parameter_estimates, parameter_bound_ranges)
Optimal_params        = [x for x in popt]
Optimal_params_err    = [x for x in perr]
Optimal_params_names  = ["$\\eta$", "Constant C"]
# Optimal_params_names  = ["$\\xi$", "Constant A"]

chi_square_min        = chi_square(function, xdata_fit, ydata_fit, yerror_fit, *Optimal_params)
Optimal_params.append(chi_square_min)
Optimal_params_err.append("")
Optimal_params_names.append("Min Chi Square")


xfit = xdata_fit
yfit = function(xfit, *Optimal_params)

title    = "Spin Correlation near T_c ($T = " + str(T_SC[iteration]) + " (J/k_B)$) vs. Distance (Power Law Fit)"
y_label  = "Spin Correlation $R(x)$ (Unitless)"
x_label  = "Distance $x$ (Spacing)"
fit_eq   = "$R(x) = C*|x|^{-\\eta}$"

fig_num = plot_fit(xdata_data, ydata_data, yerror_data, xfit, yfit, Optimal_params, Optimal_params_err, Optimal_params_names, fig_num, title=title, x_label=x_label, y_label=y_label, fit_eq=fit_eq)
fig_num = plot_residuals(xdata_fit, ydata_fit, yerror_fit, xfit, yfit, fig_num, title=title, x_label=x_label, y_label=y_label, log_scale_y=True, log_scale_x=True)
print("Spin Correlation Near T_c: END\n")
plt.show()
#-----------------------------------------------------------------------------------------------------------------------------------------



