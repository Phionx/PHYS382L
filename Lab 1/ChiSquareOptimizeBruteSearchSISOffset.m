clear all;
close all;
%Setup File Name
clear;
date_taken              = '20190221'; %Date
date_written            = '20190221'; %Date Fitted
tunneling_type          =      'SIS'; %1: SIN, 2: NIN, 3: SIS
junction_type           =        '2'; %1,2, or 3
trial                   =        '0'; %Just so we can cleanly store data, positive trials are forward sweeps, negative trials are backwards sweeps
writing                 =          0; %if 1, then save csv and images, anything else = no save

%File Names
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
%           [Delta_al (meV), T (K),  R_0 (Ohms), Delta_pb (meV), Offset] 
L         = [           .02,   1.2,  R_0_fit-5,             1.3,  ];%LOWER BOUND
U         = [           .09,   1.4,  R_0_fit+5, 1.5];%UPPER BOUND
%increment = [            .002,   .05,         .2];%INCREMENT SIZE
%NUMBER OF PARAMETER POINTS ON EACH AXIS [Delta_al, T, R_0]
% len       = zeros(1,3);
% for k=1:4
%     len(k) = round(round((U(k)-L(k))/increment(k))+1);
% end

%len       = [             21,     21,          21, 21];
len        = [             11,     5,         5, 11];

%Actual Parameter Space 
parameter_grid_delta_al = linspace(L(1,1), U(1,1),len(1)); 
parameter_grid_T        = linspace(L(1,2), U(1,2),len(2)); 
parameter_grid_R_0      = linspace(L(1,3), U(1,3),len(3));
parameter_grid_delta_pb = linspace(L(1,4), U(1,4),len(4));

%Chisquare Values over Parameter Space
chisquare_grid       = zeros(len(1), len(2), len(3), len(4));

disp(len);
for k = 1:len(1)%Vary delta_al
    disp(k);
    for j = 1:len(2)%Vary T
        for l = 1:len(3)%Vary R_0
            for m = 1:len(4)%Vary delta_pb
                chisquare_grid(k,j,l,m) = chisquare([parameter_grid_delta_al(k), parameter_grid_T(j), parameter_grid_R_0(l), parameter_grid_delta_pb(m)]);
            end
        end
    end
end


%Plotting Best Fit
%--------------------------------------------------------------------------
%Theoretical Fit
%Best Fit Paramaters
[minimum, min_index]              = min(chisquare_grid(:));
[delta_al_index, T_index, R_0_index, delta_pb_index] = ind2sub(size(chisquare_grid),min_index);
delta_al                             = parameter_grid_delta_al(delta_al_index);
T                                 = parameter_grid_T(T_index);
R_0                               = parameter_grid_R_0(R_0_index);
delta_pb                          = parameter_grid_delta_pb(delta_pb_index);
Fits                              = zeros(measurement_length, 1);%So that we can export Fit data
Fits(1:4, 1)                      = [delta_al, T, R_0, delta_pb];

%Minimum Chisquare
[chisquare_min_val, itot_min]     = SIScurr(delta_al, R_0, T, delta_pb, Input_V_j', Current_I_j', Total_Error_Current_I_j');
Fits(5,1)                         = chisquare_min_val;
%Best Fit Parameter Estimate Uncertainties
Uncertainty_Index                 = zeros(2, 4); %min_index, max_index = index in parameter space where chi_square doubles for 3 parameters
Uncertainties                     = zeros(2, 4); %min,max = 2, numParameters = 3
Uncertainty_Avg                   = zeros(measurement_length,1);%padded with xeros so that we can export it 

%Delta_al Uncertainty
Uncertainty_Index(1,1)          = 1;
Uncertainty_Index(2,1)          = len(1);


for k= fliplr(1:delta_al_index)
    if chisquare_grid(k, T_index, R_0_index, delta_pb_index)/chisquare_min_val >= 2
        Uncertainty_Index(1,1)  = k;
        break
    end
end
for k = delta_al_index:len(1)
    if chisquare_grid(k, T_index, R_0_index, delta_pb_index)/chisquare_min_val >= 2
        Uncertainty_Index(2,1)  = k;
        break
    end
end

Uncertainties(1,1)  = abs(parameter_grid_delta_al(delta_al_index) - parameter_grid_delta_al(Uncertainty_Index(1,1)));
Uncertainties(2,1)  = abs(parameter_grid_delta_al(delta_al_index) - parameter_grid_delta_al(Uncertainty_Index(2,1)));
Uncertainty_Avg(1)  = (Uncertainties(1,1) + Uncertainties(2,1))/2;

%T Uncertainty
Uncertainty_Index(1,2)          = 1;
Uncertainty_Index(2,2)          = len(2);
for k = fliplr(1:T_index)
    if chisquare_grid(delta_al_index, k, R_0_index, delta_pb_index)/chisquare_min_val >= 2
        Uncertainty_Index(1,2)  = k;
        break
    end
end
for k = T_index:len(2)
    if chisquare_grid(delta_al_index, k, R_0_index, delta_pb_index)/chisquare_min_val >= 2
        Uncertainty_Index(2,2)  = k;
        break
    end
end


Uncertainties(1,2)  = abs(parameter_grid_T(T_index) - parameter_grid_T(Uncertainty_Index(1,2)));
Uncertainties(2,2)  = abs(parameter_grid_T(T_index) - parameter_grid_T(Uncertainty_Index(2,2)));
Uncertainty_Avg(2)  = (Uncertainties(1,2) + Uncertainties(2,2))/2;

%R_0 Uncertainty
Uncertainty_Index(1,3)          = 1;
Uncertainty_Index(2,3)          = len(3);
for k = fliplr(1:R_0_index)
    if chisquare_grid(delta_al_index, T_index, k, delta_pb_index)/chisquare_min_val >= 2
        Uncertainty_Index(1,3)  = k;
        break
    end
end
for k = R_0_index:len(3)
    if chisquare_grid(delta_al_index, T_index, k, delta_pb_index)/chisquare_min_val >= 2
        Uncertainty_Index(2,3)  = k;
        break
    end
end

Uncertainties(1,3)  = abs(parameter_grid_R_0(R_0_index) - parameter_grid_R_0(Uncertainty_Index(1,3)));
Uncertainties(2,3)  = abs(parameter_grid_R_0(R_0_index) - parameter_grid_R_0(Uncertainty_Index(2,3)));
Uncertainty_Avg(3)  = (Uncertainties(1,3) + Uncertainties(2,3))/2;

%delta_pb Uncertainty
Uncertainty_Index(1,4)          = 1;
Uncertainty_Index(2,4)          = len(4);
for k = fliplr(1:delta_pb_index)
    if chisquare_grid(delta_al_index, T_index, R_0_index, k)/chisquare_min_val >= 2
        Uncertainty_Index(1,4)  = k;
        break
    end
end
for k = delta_pb_index:len(4)
    if chisquare_grid(delta_al_index, T_index, R_0_index, k)/chisquare_min_val >= 2
        Uncertainty_Index(2,4)  = k;
        break
    end
end

Uncertainties(1,4)  = abs(parameter_grid_delta_pb(delta_pb_index) - parameter_grid_delta_pb(Uncertainty_Index(1,4)));
Uncertainties(2,4)  = abs(parameter_grid_delta_pb(delta_pb_index) - parameter_grid_delta_pb(Uncertainty_Index(2,4)));
Uncertainty_Avg(4)  = (Uncertainties(1,4) + Uncertainties(2,4))/2;

%Plotting and Fitting
%--------------------------------------------------------------------------
%Linear Fit to find Approximate Resistance
%trial_length = 255;%Number of Points in each Trial





%FIGURE 1 (MAIN)
x        = Input_V_j;
y        = Current_I_j;
error    = Total_Error_Current_I_j;
theory_y = itot_min;
figure(1);
%Linear Fit for R_0
plot(x, Linear_Fit(1)*x + Linear_Fit(2), 'b-.')
hold on;
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
equation = strcat(equation, '\newline \Delta_{al}: ', string(delta_al),  ' ', char(177), ' ', string(Uncertainty_Avg(1)), ' meV');
equation = strcat(equation, '\newline T: ', string(T),      ' ', char(177), ' ', string(Uncertainty_Avg(2)), ' K');
equation  = strcat(equation, '\newline R_0: ', string(R_0),      ' ', char(177), ' ', string(Uncertainty_Avg(3)), ' \Omega');
equation  = strcat(equation, '\newline \Delta_{pb}: ', string(delta_pb),      ' ', char(177), ' ', string(Uncertainty_Avg(4)), ' \Omega');
yL=get(gca,'YLim'); 
xL=get(gca,'XLim');   
%text((xL(1)+xL(2))/2,(yL(1) + yL(2))*1/2, equation,...
text(0.0015,-0.00001,equation,...
      'HorizontalAlignment','left',...
      'VerticalAlignment','top',...
      'BackgroundColor',[.8 .8 .8],...
      'FontSize',12);

 
disp(chisquare_min_val);

disp(delta_al);
disp(T);
disp(R_0);
disp(delta_pb);

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
plot(x, Linear_Fit(1)*x + Linear_Fit(2), 'b-.')
hold on;
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
    delta_al = x(1,1);%meV
    T        = x(1,2);%K
    R_0      = x(1,3);%Ohms
    delta_pb = x(1,4);%meV
    Voltages = Input_V_j;
    Currents = Current_I_j;
    Error_Currents = Total_Error_Current_I_j;
    [chisquare_val, itot] = SIScurr(delta_al, R_0, T, delta_pb, Voltages', Currents', Error_Currents');
end

