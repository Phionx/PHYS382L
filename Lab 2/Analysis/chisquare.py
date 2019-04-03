import numpy as np
import scipy.optimize as opt
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import math
import csv


#DATA PARSING
#-----------------------------------------------------------------------------------------------------------------------------------------
#T,E,M,C,X  - Temperature, Energy, Magnetism, Specific Heat, Susceptibility
EMfile = 'Data/APR:1:2019/4_01_N100_1.20T3.00_EM_v0.csv'
with open(EMfile, 'r') as csvEMFile:
	readerEM = csv.reader(csvEMFile)
	readerEM   = list(readerEM)
csvEMFile.close()
readerEM  = readerEM[4:]
# first_col = [row[0] for row in readerEM]
# max_ind   = first_col.index('')
# readerEM  = readerEM[:max_ind]

#Data Parsing
T_EM   = np.array([row[0] for row in readerEM]).astype(np.float)
E_mean = np.array([row[2] for row in readerEM]).astype(np.float)
E_std  = np.array([row[3] for row in readerEM]).astype(np.float)
M_mean = np.array([row[4] for row in readerEM]).astype(np.float)
M_mean = np.abs(M_mean)
M_std  = np.array([row[5] for row in readerEM]).astype(np.float)
C_mean = np.array([row[6] for row in readerEM]).astype(np.float)
C_std  = np.array([row[7] for row in readerEM]).astype(np.float)
X_mean = np.array([row[8] for row in readerEM]).astype(np.float)
X_std  = np.array([row[9] for row in readerEM]).astype(np.float)

#Spin Correlations
SCfile = 'Data/APR:1:2019/4_01_N100_1.20T3.00_SC_v0.csv'
with open(SCfile, 'r') as csvSCFile:
	readerSC = csv.reader(csvSCFile)
	readerSC = list(readerSC)
csvSCFile.close()
d_max     = int(readerSC[1][0])/2-1 #N value/2
readerSC  = readerSC[4:]
# first_col = [row[0] for row in readerSC]
# readerSC  = readerSC[:max_ind]

#Data Parsing
T_SC    = np.array([row[0] for row in readerSC]).astype(np.float)
mean_SC = np.array([row[2:d_max+2] for row in readerSC]).astype(np.float)
std_SC  = np.array([row[d_max+2:2*d_max+2] for row in readerSC]).astype(np.float)
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
	scale = parameters[0]
	add   = parameters[1]

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
	return scale*ans + add

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
	plt.plot(xdata, ydata, 'k.', color='#3F7F4C', label="Simulation Data")
	plt.fill_between(xdata, ydata-yerror, ydata+yerror, alpha=.8, edgecolor='#3F7F4C', facecolor='#7EFF99', linewidth=0)

	#Best Theory
	theory_label = "Theory"

	
	theory_label += " " + fit_eq

	theory_label += "\n" + params_names[-1] + ": " + str(params[-1])

	for i in range(len(params)-1):
		theory_label += "\n"
		theory_label += params_names[i] + ": "
		theory_label += str(params[i]) + " pm "
		theory_label += str(params_err[i]) 

	print(xtheory)
	print(ytheory)
	plt.plot(xtheory, ytheory, 'k-', color='#CC4F1B', label=theory_label)
	plt.legend(loc='upper right')
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	plt.title(title)

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

	# print("residuals(data): " + str(ydata))
	# print("residuals(theory): " + str(ytheory))
	yres    = [(ydata[i] - ytheory[i])/yerror[i] for i in range(data_points)]
	xres    = xdata
	plt.figure(fig_num)
	plt.plot(xres, yres, 'k-', color='#3F7F4C')
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	plt.title(title)
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
scale_constant            = 1.0
add_constant              = 1.0
parameter_estimates       = [scale_constant, add_constant]

var                       = [[-100, -100],[100, 100]]
parameter_bound_ranges    = ([parameter_estimates[i] + var[0][i] for i in range(len(var[0]))], [parameter_estimates[i] + var[1][i] for i in range(len(var[1]))])
xdata_data                = T_EM
ydata_data                = E_mean
yerror_data               = E_std
function                  = energy_func
fig_num                   = 1

xdata_fit                 = xdata_data
ydata_fit                 = ydata_data
yerror_fit                = yerror_data


popt, perr            = optimal_fit(function, xdata_fit, ydata_fit, yerror_fit, parameter_estimates, parameter_bound_ranges)
Optimal_params        = [x for x in popt]
Optimal_params_err    = [x for x in perr]
Optimal_params_names  = ["Scale Constant A", "Add Constant B"]

chi_square_min        = chi_square(function, xdata_fit, ydata_fit, yerror_fit, *Optimal_params)
Optimal_params.append(chi_square_min)
Optimal_params_names.append("Min Chi Square")


xfit = xdata_fit
yfit = function(xfit, *Optimal_params)

title    = "Energy vs. Temperature"
y_label  = "Energy (Units)"
x_label  = "Temperature (K)"
fit_eq   = "$E = AE_\\{original\\}(t) + B$"

plot_fit(xdata_data, ydata_data, yerror_data, xfit, yfit, Optimal_params, Optimal_params_err, Optimal_params_names, fig_num, title=title, x_label=x_label, y_label=y_label, fit_eq=fit_eq)

fig_num += 1
title    = "Energy vs. Temperature (Residuals)"
y_label  = "Energy (Units)"
x_label  = "Temperature (K)"
plot_residuals(xdata_fit, ydata_fit, yerror_fit, xfit, yfit, fig_num, title=title, x_label=x_label, y_label=y_label)
#-----------------------------------------------------------------------------------------------------------------------------------------

#Magnetization Full Fit
#-----------------------------------------------------------------------------------------------------------------------------------------
Estimate_T_c_mag          = 2.269

function                  =	mag_func
xdata_data                = T_EM
ydata_data                = M_mean
yerror_data               = M_std
fig_num                   += 1

fit_fraction_left_outer   = 1.0
fit_fraction_left_inner   = 1.0/30

fit_fraction_right_outer  = 1.0
fit_fraction_right_inner  = 1.0/58

fit_center_index          = find_nearest(xdata_data, Estimate_T_c_mag)
data_points               = len(xdata_data)
xdata_fit_left            = xdata_data[fit_center_index - int(fit_fraction_left_outer*fit_center_index):fit_center_index - int(fit_fraction_left_inner*fit_center_index)]
xdata_fit_right           = xdata_data[fit_center_index + int(fit_fraction_right_inner*(data_points-fit_center_index-1)):fit_center_index + int(fit_fraction_right_outer*(data_points-fit_center_index))]
xdata_fit                 = np.concatenate((xdata_fit_left, xdata_fit_right), axis=None)
ydata_fit_left            = ydata_data[fit_center_index - int(fit_fraction_left_outer*fit_center_index):fit_center_index - int(fit_fraction_left_inner*fit_center_index)]
ydata_fit_right           = ydata_data[fit_center_index + int(fit_fraction_right_inner*(data_points-fit_center_index-1)):fit_center_index + int(fit_fraction_right_outer*(data_points-fit_center_index))]
ydata_fit                 = np.concatenate((ydata_fit_left, ydata_fit_right), axis=None)

yerror_fit_left            = yerror_data[fit_center_index - int(fit_fraction_left_outer*fit_center_index):fit_center_index - int(fit_fraction_left_inner*fit_center_index)]
yerror_fit_right           = yerror_data[fit_center_index + int(fit_fraction_right_inner*(data_points-fit_center_index-1)):fit_center_index + int(fit_fraction_right_outer*(data_points-fit_center_index))]
yerror_fit                 = np.concatenate((yerror_fit_left, yerror_fit_right), axis=None)


Optimal_params        = []
Optimal_params_err    = []
Optimal_params_names  = []

chi_square_min        = chi_square(function, xdata_fit, ydata_fit, yerror_fit)
Optimal_params.append(chi_square_min)
Optimal_params_names.append("Min Chi Square")


xfit = xdata_fit
yfit = [function(x) for x in xfit]


title    = "Magnetization vs. Temperature"
y_label  = "Magnetization (Units)"
x_label  = "Temperature (K)"
fit_eq   = ""

plot_fit(xdata_data, ydata_data, yerror_data, xfit, yfit, Optimal_params, Optimal_params_err, Optimal_params_names, fig_num, title=title, x_label=x_label, y_label=y_label, fit_eq=fit_eq)

fig_num += 1
title    = "Magnetization vs. Temperature (Residuals)"
y_label  = "Magnetization (Units)"
x_label  = "Temperature (K)"
plot_residuals(xdata_fit, ydata_fit, yerror_fit, xfit, yfit, fig_num, title=title, x_label=x_label, y_label=y_label)
#-----------------------------------------------------------------------------------------------------------------------------------------


#Magnetization: T_C, \beta
#-----------------------------------------------------------------------------------------------------------------------------------------
Estimate_T_c_beta         = 2.269
Estimate_beta             = 1.0/8
Constant                  = 1.0
parameter_estimates       = [Estimate_T_c_beta, Estimate_beta, Constant]
var                       = [[-.1, 0, -1000],[.1, .2, 10000]]

parameter_bound_ranges    = ([parameter_estimates[i] + var[0][i] for i in range(len(var[0]))], [parameter_estimates[i] + var[1][i] for i in range(len(var[1]))])
function                  = beta_func
xdata_data                = T_EM
ydata_data                = M_mean
yerror_data               = M_std
fig_num                   += 1

parameter_estimates[0]    = 2.3
fit_center_index          = find_nearest(xdata_data, parameter_estimates[0])

fit_fraction_left_outer   = 2.0/4
fit_fraction_left_inner   = 1.0/10
fit_fraction_right_outer  = 0
fit_fraction_right_inner  = 0


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
Optimal_params_names  = ["$T_c$", "$\\beta$", "Constant C"]

chi_square_min        = chi_square(function, xdata_fit, ydata_fit, yerror_fit, *Optimal_params)
Optimal_params.append(chi_square_min)
Optimal_params_names.append("Min Chi Square")


xfit = xdata_fit
yfit = function(xfit, *Optimal_params)

title    = "Magnetization vs. Temperature"
y_label  = "Magnetization (Units)"
x_label  = "Temperature (K)"
fit_eq   = "$|M| = C|t|^{\\beta}$"

plot_fit(xdata_data, ydata_data, yerror_data, xfit, yfit, Optimal_params, Optimal_params_err, Optimal_params_names, fig_num, title=title, x_label=x_label, y_label=y_label, fit_eq=fit_eq)

fig_num += 1
title    = "Magnetization vs. Temperature (Residuals)"
y_label  = "Magnetization (Units)"
x_label  = "Temperature (K)"
plot_residuals(xdata_fit, ydata_fit, yerror_fit, xfit, yfit, fig_num, title=title, x_label=x_label, y_label=y_label)
#-----------------------------------------------------------------------------------------------------------------------------------------



#Susceptibility: T_C, \gamma
#-----------------------------------------------------------------------------------------------------------------------------------------
Estimate_T_c_gamma        = 2.269
Estimate_gamma            = 7.0/4
Constant                  = 1.0
parameter_estimates       = [Estimate_T_c_gamma, Estimate_gamma, Constant]
var                       = [[-.01, 0, -1000],[.1, .01, 10000]]

parameter_bound_ranges    = ([parameter_estimates[i] + var[0][i] for i in range(len(var[0]))], [parameter_estimates[i] + var[1][i] for i in range(len(var[1]))])
function                  = gamma_func
xdata_data                = T_EM
ydata_data                = X_mean
yerror_data               = X_std
fig_num                   += 1



# parameter_estimates[0]    = 2.3
fit_center_index          = find_nearest(xdata_data, parameter_estimates[0])

fit_fraction_left_outer   = 0
fit_fraction_left_inner   = 0
fit_fraction_right_outer  = 1.0/3
fit_fraction_right_inner  = 1.0/12

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
Optimal_params_names  = ["$T_c$", "$\\gamma$", "Constant C"]

chi_square_min        = chi_square(function, xdata_fit, ydata_fit, yerror_fit, *Optimal_params)
Optimal_params.append(chi_square_min)
Optimal_params_names.append("Min Chi Square")


xfit = xdata_fit
yfit = function(xfit, *Optimal_params)

title    = "Susceptibility vs. Temperature"
y_label  = "Susceptibility (Units)"
x_label  = "Temperature (K)"
fit_eq   = "$\\chi = C|t|^{-\\gamma}$"

plot_fit(xdata_data, ydata_data, yerror_data, xfit, yfit, Optimal_params, Optimal_params_err, Optimal_params_names, fig_num, title=title, x_label=x_label, y_label=y_label, fit_eq=fit_eq)

fig_num += 1
title    = "Susceptibility vs. Temperature (Residuals)"
y_label  = "Susceptibility (Units)"
x_label  = "Temperature (K)"
plot_residuals(xdata_fit, ydata_fit, yerror_fit, xfit, yfit, fig_num, title=title, x_label=x_label, y_label=y_label)
#-----------------------------------------------------------------------------------------------------------------------------------------

#Specific Heat: T_c, \alpha
#-----------------------------------------------------------------------------------------------------------------------------------------
Estimate_T_c_alpha        = 2.269
Estimate_alpha            = 7.0/4
Constant                  = 1.0
parameter_estimates       = [Estimate_T_c_alpha, Estimate_alpha, Constant]
var                       = [[-.01, 0, -1000],[.01,.01, 10000]]

parameter_bound_ranges    = ([parameter_estimates[i] + var[0][i] for i in range(len(var[0]))], [parameter_estimates[i] + var[1][i] for i in range(len(var[1]))])
function                  = alpha_func
xdata_data                = T_EM
ydata_data                = C_mean
yerror_data               = C_std
fig_num                   += 1



# parameter_estimates[0]    = 2.3
fit_center_index          = find_nearest(xdata_data, parameter_estimates[0])

fit_fraction_left_outer   = 1.0/6
fit_fraction_left_inner   = 1.0/10
fit_fraction_right_outer  = 0.0/3
fit_fraction_right_inner  = 0.0/10



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
Optimal_params_names  = ["$T_c$", "$\\alpha$", "Constant C"]

chi_square_min        = chi_square(function, xdata_fit, ydata_fit, yerror_fit, *Optimal_params)
Optimal_params.append(chi_square_min)
Optimal_params_names.append("Min Chi Square")


xfit = xdata_fit
yfit = function(xfit, *Optimal_params)

title    = "Specific Heat vs. Temperature"
y_label  = "Specific Heat (Units)"
x_label  = "Temperature (K)"
fit_eq   = "$\\chi = C|t|^{-\\alpha}$"

plot_fit(xdata_data, ydata_data, yerror_data, xfit, yfit, Optimal_params, Optimal_params_err, Optimal_params_names, fig_num, title=title, x_label=x_label, y_label=y_label, fit_eq=fit_eq)

fig_num += 1
title    = "Specific Heat vs. Temperature (Residuals)"
y_label  = "Specific Heat (Units)"
x_label  = "Temperature (K)"
plot_residuals(xdata_fit, ydata_fit, yerror_fit, xfit, yfit, fig_num, title=title, x_label=x_label, y_label=y_label)

plt.show()
#-----------------------------------------------------------------------------------------------------------------------------------------
