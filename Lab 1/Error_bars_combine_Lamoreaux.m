%This matlab file is to combine the error bars of our junction current I_j
%and our junction voltage V_j. To do this, we must approximately linearly
%fit each of our curves and propogate error in V_j using the relation 
% I_j = m*V_j, where m = 1/R_effective, where R_effective is the
% approximate linear resistance of our junction


%Setup File Name
date_acquired      = '20190219'; %Date of Data Acquisition
date_written       = '20190219'; %Date of Data Analysis
tunneling_type     = 'SIS';%1: SIN, 2: NIN, 3: SIS
junction_type      = '3';%1,2, or 3
trial              = '2';%Made so we can cleanly store data, positive trials are forward sweeps, negative trials are backwards sweeps


file_reading              = strcat('measurementsAnalysis/LamoreauxErrorBars/', date_acquired, '_', junction_type, '_', tunneling_type, '_', 'Trial', string(trial), '.csv');
data                      = csvread(file_reading);
Current_I_j               = data(:, 1)';
Error_Current_I_j         = data(:, 2)';
Input_V_j                 = data(:, 3)';
Error_Input_V_j           = data(:,4)';

file_writing              = strcat('measurementsAnalysis/LamoreauxCombinedErrorBars/', date_written , '_', junction_type, '_', tunneling_type, '_', 'Trial', string(trial), '.csv');

%Plotting and Fitting
%--------------------------------------------------------------------------
%Linear Fit to find Approximate Resistance
scatter(Input_V_j, Current_I_j)
Fit = polyfit(Input_V_j, Current_I_j,1);
hold on;

xlabel('Junction Voltage (V)');
ylabel('Junction Current (A)');
title('I(V) Curve');

%Linear Fit Equation
equation = sprintf('I = %.6f V + %.6f', Fit(1), Fit(2));
%text(.001, 0, equation, 'FontSize', 10);
yL=get(gca,'YLim'); 
xL=get(gca,'XLim');   
text((xL(1)+xL(2))/2,yL(2),equation,...
      'HorizontalAlignment','left',...
      'VerticalAlignment','top',...
      'BackgroundColor',[1 1 1],...
      'FontSize',12);
plot(Input_V_j, Fit(1)*Input_V_j + Fit(2), 'r-.');

%Combine Horizontal and Veritical Error Bars
%--------------------------------------------------------------------------
measurement_length = length(Error_Current_I_j);
%R_effective   = 1/Fit(1); %Approximate I = (1/R)V, R_eff = 1/m
R_effective    = 1/.009934;
Converted_Error_Voltage_V_j = (1/R_effective)*Error_Input_V_j;
Total_Error_Current_I_j = zeros(1,measurement_length);
i = 1;
while i <= measurement_length
    Total_Error_Current_I_j(i)     = sqrt(Error_Current_I_j(i)^2 + Converted_Error_Voltage_V_j(i)^2);
    i = i + 1;
end

%Data Export
%--------------------------------------------------------------------------
export_data = [Input_V_j; Current_I_j; Total_Error_Current_I_j];
export_data = export_data';

%Export to Right File
dlmwrite(file_writing, export_data, 'delimiter', ',', '-append');

