clear all;
close all;
%Setup File Name
clear;
date_taken              = '20190220'; %Date
date_written            = '20190220'; %Date Fitted
tunneling_type          =      'SIN'; %1: SIN, 2: NIN, 3: SIS
junction_type           =        '2'; %1,2, or 3
trial                   =        '0'; %Just so we can cleanly store data, positive trials are forward sweeps, negative trials are backwards sweeps

%File Names
file_read               = strcat('Data/measurementsAnalysis/FittedData/', date_taken, '_', junction_type, '_', tunneling_type, '_', 'Trial', string(trial), '.csv');
file_write_image        = strcat('Figures/Residuals/', date_written, '_', junction_type, '_', tunneling_type);
data                    = csvread(file_read);
global Input_V_j;
global Current_I_j;
global Total_Error_Current_I_j;
global Theory_Current_I_j;
Input_V_j               = data(:, 1);
Current_I_j             = data(:, 2);
Total_Error_Current_I_j = data(:, 3);
Theory_Current_I_j      = data(:, 4);
measurement_length      = length(Input_V_j);


%Calculate Residuals
Residuals_fit           = zeros(measurement_length,1);
for k=1:measurement_length
    Residuals_fit(k)        = (Current_I_j(k) - Theory_Current_I_j(k))/(Total_Error_Current_I_j(k));
end

%FIGURE 1 (RESIDUALS)
x        = Input_V_j;
y        = Residuals_fit;
%Residuals
plot(x, y, 'g')
hold on;


axis([-inf inf -inf inf]);
%axis([0 inf 0 inf]);
xlabel('Junction Voltage (V)');
ylabel('Residuals');
title(strcat(tunneling_type, ' Residual Curve, Junction ', junction_type));

print(file_write_image, '-dpng');

