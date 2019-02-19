function [temp] = resistance_thermometer(lockin_output)

LHe_res = 1425;
LHe_temp = 4.2;
LN2_res = 1047;
LN2_temp = 77.4;
therm_slope = log10(LHe_temp/LN2_temp)/log10(LHe_res/LN2_res);

% values used to calibrate resistance thermometer (Ben Brubaker, 4/18/2014),
%       assuming log(R) vs. log(T) is linear between 1 K and 100 K

temp_circuit_r = 131000; % fixed current source resistance for thermometry four-wire measurement
excitation_voltage = 0.01; % current source excitation voltage

therm_r = lockin_output*(temp_circuit_r/excitation_voltage);

temp = LHe_temp*(therm_r/LHe_res)^therm_slope;