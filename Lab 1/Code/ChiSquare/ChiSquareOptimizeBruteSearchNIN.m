clear all;
close all;
%Setup File Name
clear;
date_taken              = '20190221'; %Date
date_written            = '20190221'; %Date Fitted
tunneling_type          =      'NIN'; %1: SIN, 2: NIN, 3: SIS
junction_type           =        '1'; %1,2, or 3
trial                   =        '1'; %Just so we can cleanly store data, positive trials are forward sweeps, negative trials are backwards sweeps
writing                 =          0; %if 1, then save csv and images, anything else = no save

%File Names
file_read               = strcat('Data/measurementsAnalysis/SeparateNIN/', date_taken, '_', junction_type, '_', tunneling_type, '_', 'Trial', string(trial), '.csv');
%file_read               = strcat('Data/measurementsAnalysis/SortedWrongAnalysis/', date_taken, '_', junction_type, '_', tunneling_type, '_', 'Trial', string(trial), '.csv');
%file_read               = strcat('Data/measurementsAnalysis/FinalData/', date_taken, '_', junction_type, '_', tunneling_type, '_', 'Trial', string(trial), '.csv');
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


%Minimization
%--------------------------------------------------------------------------
%Linear Fit on End Points to Find R_0
disp(Input_V_j);
end_points        = zeros(2,7);
end_points(1,1:3)  = Input_V_j(1:3,1);
end_points(1,4:6)  = Input_V_j((measurement_length - 2):measurement_length,1);
end_points(2,1:3)  = Current_I_j(1:3,1);
end_points(2,4:6)  = Current_I_j((measurement_length - 2):measurement_length,1);
end_points(1,7)   = 0;
end_points(2,7)   = 0;

%Linear Fit on End Points
Linear_Fit = polyfit(end_points(1,:), end_points(2,:), 1);
R_0_fit = 1/Linear_Fit(1);  



%Parameter Ranges
%           [R_0] 
L         = [R_0_fit-3];%LOWER BOUND
U         = [R_0_fit+3];%UPPER BOUND
len       = [     1000];
%Actual Parameter Space 
parameter_grid_R_0   = linspace(L(1,1), U(1,1),len(1));

%Chisquare Values over Parameter Space
chisquare_grid       = zeros(len(1),1);

for l = 1:len(1)%Vary Delta
    disp(l);
    chisquare_grid(l) = chisquare([parameter_grid_R_0(l)]);
end


%Plotting Best Fit
%--------------------------------------------------------------------------
%Theoretical Fit
%Best Fit Paramaters
[minimum, min_index]              = min(chisquare_grid(:));
[R_0_index]                       = ind2sub(size(chisquare_grid),min_index);
R_0                               = parameter_grid_R_0(R_0_index);
Fits                              = zeros(measurement_length, 1);%So that we can export Fit data
Fits(1, 1)                        = R_0;

%Minimum Chisquare
[chisquare_min_val, itot_min]     = NINcurr(R_0, Input_V_j', Current_I_j', Total_Error_Current_I_j');
Fits(2,1)                         = chisquare_min_val;
%Best Fit Parameter Estimate Uncertainties
Uncertainty_Index                 = zeros(2, 1); %min_index, max_index = index in parameter space where chi_square doubles for 3 parameters
Uncertainties                     = zeros(2, 1); %min,max = 2, numParameters = 3
Uncertainty_Avg                   = zeros(measurement_length,1);%padded with xeros so that we can export it 


%R_0 Uncertainty
Uncertainty_Index(1,1)          = 1;
Uncertainty_Index(2,1)          = len(1);
for k = fliplr(1:R_0_index)
    if chisquare_grid(k)/chisquare_min_val >= 2
        Uncertainty_Index(1,1)  = k;
        break
    end
end
for k = R_0_index:len(1)
    if chisquare_grid(k)/chisquare_min_val >= 2
        Uncertainty_Index(2,1)  = k;
        break
    end
end

Uncertainties(1,1)  = abs(parameter_grid_R_0(R_0_index) - parameter_grid_R_0(Uncertainty_Index(1,1)));
Uncertainties(2,1)  = abs(parameter_grid_R_0(R_0_index) - parameter_grid_R_0(Uncertainty_Index(2,1)));
Uncertainty_Avg(1)  = (Uncertainties(1,1) + Uncertainties(2,1))/2;


%Plotting and Fitting
%--------------------------------------------------------------------------
%FIGURE 1 (MAIN)
x        = Input_V_j;
y        = Current_I_j;
error    = Total_Error_Current_I_j;
theory_y = itot_min;
figure(1);
%Linear Fit for R_0
%plot(x, Linear_Fit(1)*x + Linear_Fit(2), 'b-.')
%hold on;
%Data
shadedErrorBar(x, y, error, 'lineprops', '-r','transparent',false,'patchSaturation',0.075);
hold on;
%Theory Fit
plot(x, theory_y, 'g')


axis([-inf inf -inf inf]);
%axis([0 inf 0 inf]);
xlabel('Junction Voltage (V)');
ylabel('Junction Current (A)');
title(strcat(tunneling_type, ' I(V) Curve, Junction ', junction_type));


equation = strcat('\chi^2_{\nu}: ', string(chisquare_min_val));
equation  = strcat(equation, '\newline R_0: ', string(R_0),      ' ', char(177), ' ', string(Uncertainty_Avg(1)), ' \Omega');
yL=get(gca,'YLim'); 
xL=get(gca,'XLim');   
%text((xL(1)+xL(2))/2,(yL(1) + yL(2))*1/2, equation,...
text(0.002,-0.00002,equation,...
      'HorizontalAlignment','left',...
      'VerticalAlignment','top',...
      'BackgroundColor',[.8 .8 .8],...
      'FontSize',12);

 
disp(chisquare_min_val);
disp(R_0);

if writing == 1
    print(file_write_image, '-dpng');
end

%FIGURE 1 (INSET)
start_inset = measurement_length*20/40;
end_inset   = measurement_length*21/40;
start_inset = round(start_inset);
end_inset   = round(end_inset);


x        = Input_V_j(start_inset:end_inset);
y        = Current_I_j(start_inset:end_inset);
error    = Total_Error_Current_I_j(start_inset:end_inset);
theory_y = itot_min(start_inset:end_inset);
figure(2);
%Linear Fit for R_0
% plot(x, Linear_Fit(1)*x + Linear_Fit(2), 'b-.')
% hold on;
%Data
shadedErrorBar(x, y, error, 'lineprops', '-r','transparent',false,'patchSaturation',0.075);
hold on;
%Theory Fit
plot(x, theory_y, 'g')
hold on;


axis([-inf inf -inf inf]);
%axis([0 inf 0 inf]);
xlabel('Junction Voltage (V)');
ylabel('Junction Current (A)');
title(strcat(tunneling_type, ' I(V) Curve, Junction ', junction_type));



if writing == 1
    print(file_write_image_inset, '-dpng');
end

%Data Export
%--------------------------------------------------------------------------
itot_min    = itot_min';
export_data = [Input_V_j'; Current_I_j'; Total_Error_Current_I_j'; itot_min'; Fits'; Uncertainty_Avg'];
export_data = export_data';

%Export to Right File
if writing == 1
    dlmwrite(file_write_fit, export_data, 'delimiter', ',');
end

function [chisquare_val] = chisquare(x)
    global Input_V_j;
    global Current_I_j;
    global Total_Error_Current_I_j;
    R_0      = x(1,1);%Ohms
    Voltages = Input_V_j;
    Currents = Current_I_j;
    Error_Currents = Total_Error_Current_I_j;
    [chisquare_val, itot] = NINcurr(R_0, Voltages', Currents', Error_Currents');
end

