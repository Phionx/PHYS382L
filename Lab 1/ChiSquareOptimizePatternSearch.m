clear all;
close all;
%Setup File Name
clear;
date_taken              = '20190221'; %Date
date_written            = '20190221'; %Date Fitted
tunneling_type          =      'SIS'; %1: SIN, 2: NIN, 3: SIS
junction_type           =        '3'; %1,2, or 3
trial                   =        '0'; %Just so we can cleanly store data, positive trials are forward sweeps, negative trials are backwards sweeps
writing                 =          0; %if 1, then save csv and images, anything else = no save

%File Names
%file_read               = strcat('Data/measurementsAnalysis/SortedWrongAnalysis/', date_taken, '_', junction_type, '_', tunneling_type, '_', 'Trial', string(trial), '.csv');
file_read               = strcat('Data/measurementsAnalysis/FinalDataSISOffset/', date_taken, '_', junction_type, '_', tunneling_type, '_', 'Trial', string(trial), '.csv');
file_write_image        = strcat('Figures/FitFigures/', date_written, '_', junction_type, '_', tunneling_type);
file_write_image_inset  = strcat('Figures/FitFigures/', date_written, '_', junction_type, '_', tunneling_type, '_Inset');
file_write_fit          = strcat('Data/measurementsAnalysis/FittedData/', date_written, '_', junction_type, '_', tunneling_type, '_', 'Trial', string(trial), '.csv');
data                    = csvread(file_read);
global Input_V_j;
global Current_I_j;
global Total_Error_Current_I_j;
Input_V_j               = data(:, 1);
Current_I_j             = data(:, 2);
Total_Error_Current_I_j = data(:, 3);
measurement_length      = length(Input_V_j);





%Optimization
%--------------------------------------------------------------------------

%Initialize
fun = @chisquare;
x0 = [0.07, 1.35, 82.2, 1.4];
L = [0, 1, 70, 1.0];
U = [.3, 2, 90, 2.0];
[x, fval] = patternsearch(fun,x0,[],[],[],[],L,U);

disp(x);
disp(fval);

%Parameter Ranges
%           [Delta_al (meV), T (K),  R_0 (Ohms), Delta_pb (meV)] 
L         = [           .07,   1.35,  82.227, 1.4];%LOWER BOUND
U         = [           .09,   1.35,  82.227, 1.6];%UPPER BOUND

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
    delta_al = x(1,1);%meV
    T        = x(1,2);%K
    R_0      = x(1,3);%Ohms
    delta_pb = x(1,4);%meV
    Voltages = Input_V_j;
    Currents = Current_I_j;
    Error_Currents = Total_Error_Current_I_j;
    [chisquare_val, itot] = SIScurr(delta_al, R_0, T, delta_pb, Voltages', Currents', Error_Currents');
end