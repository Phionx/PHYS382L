%Setup File Name
clear;
date               = '20190212'; %Date
tunneling_type     = 'SIN';%1: SIN, 2: NIN, 3: SIS
junction_type      = '2';%1,2, or 3
trial              = '0';%Just so we can cleanly store data, positive trials are forward sweeps, negative trials are backwards sweeps

file               = strcat('measurements/CombinedErrorBars/', date, '_', junction_type, '_', tunneling_type, '_', 'Trial', string(trial), '.csv');
data                    = csvread(file);
global Input_V_j;
global Current_I_j;
global Total_Error_Current_I_j;
Input_V_j               = data(:, 1);
Current_I_j             = data(:, 2);
Total_Error_Current_I_j = data(:, 3);





%Optimization
%--------------------------------------------------------------------------

%Initialize
fun = @chisquare;
x0 = [0, 0];
L = [0, 0];
U = [100000000000, 100000000000];
[x, fval] = patternsearch(fun,x0,[],[],[],[],L,U);

% %Plotting and Fitting
% %--------------------------------------------------------------------------
% %Linear Fit to find Approximate Resistance
% trial_length = 255;%Number of Points in each Trial
% i  = length(Input_V_j)/trial_length;%Number of Trials in each File
% j  = 1;
% colors = ['r','b','g','y'];
% while j <= i
%     x  = Input_V_j((j-1)*trial_length+1:j*trial_length,1);
%     y  = Current_I_j((j-1)*trial_length+1:j*trial_length,1);
%     dy = Total_Error_Current_I_j((j-1)*trial_length+1:j*trial_length,1);
%     patch([x;flipud(x)],[y-dy;flipud(y+dy)],colors(j))
%     hold on;
%     j  = j + 1;
% end
% alpha(.2);
% x  = Input_V_j;
% y  = Current_I_j;
% scatter(x,y, '.')
% 
% 
% xlabel('Junction Voltage (V)');
% ylabel('Junction Current (A)');
% title(strcat(tunneling_type, ' I(V) Curve, Junction ', junction_type));


function [chisquare_val] = chisquare(x)
    global Input_V_j;
    global Current_I_j;
    global Total_Error_Current_I_j;
    delta = x(1,1);
    r_par = x(1,2);
    R_0 = 53;%Ohms, maybe change
    T = 4.2;%K
    Voltages = Input_V_j;
    Currents = Current_I_j;
    Error_Currents = Total_Error_Current_I_j;
    chisquare_val = SINcurr(delta, r_par, R_0, T, Voltages', Currents', Error_Currents');
end