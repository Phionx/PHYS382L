clear all;
close all;
%Setup File Name
clear;
date                    = '20190212'; %Date
date_written            = '20190218'; %Date Fitted
tunneling_type          =      'SIN'; %1: SIN, 2: NIN, 3: SIS
junction_type           =        '2'; %1,2, or 3
trial                   =        '0'; %Just so we can cleanly store data, positive trials are forward sweeps, negative trials are backwards sweeps

%File Names
file_read               = strcat('measurementsAnalysis/CombinedErrorBars/', date, '_', junction_type, '_', tunneling_type, '_', 'Trial', string(trial), '.csv');
file_write              = strcat('FitFigures/', date_written, '_', junction_type, '_', tunneling_type, '_', 'Trial', string(trial));

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
%Calculating R_0 Initial Estimate with Linear Fit Equation
%Finding End Points
end_points_num = measurement_length/255;%255 is num of data points in one trial, since trials are concat together, we need this to get the end points for our R_0 fit
end_points    = zeros(2,end_points_num*3);
for k=0:end_points_num-1
    end_points(1,1 + k) = Input_V_j(1+255*k);
    end_points(2,1 + k) = Current_I_j(1+255*k);
    end_points(1,2 + k) = Input_V_j(128+255*k);
    end_points(2,2 + k) = Current_I_j(128+255*k);
    end_points(1,3 + k) = Input_V_j(255+255*k);
    end_points(2,3 + k) = Current_I_j(255+255*k);
end
    
%Linear Fit on End Points
Fit = polyfit(end_points(1,:), end_points(2,:),1);
R_0_fit = 1/Fit(1);
equation = sprintf('R_0 = %.6f', R_0_fit);
%text(.001, 0, equation, 'FontSize', 10);
yL=get(gca,'YLim'); 
xL=get(gca,'XLim');   
%text((xL(1)+xL(2))/3,yL(2)*4/5,equation,...
text(0.02,0.0003,equation,...
      'HorizontalAlignment','left',...
      'VerticalAlignment','top',...
      'BackgroundColor',[1 1 1],...
      'FontSize',12);
plot(Input_V_j, Fit(1)*Input_V_j + Fit(2), 'b-.')
hold on;


%Parameter Ranges
%           [Delta_pb (meV), T (K), R_0 (Ohms)] 
L         = [            .9,   4.0,  R_0_fit-10];%LOWER BOUND
U         = [           1.3,   4.6,  R_0_fit+10];%UPPER BOUND
increment = [           .02,   .05,         .2];%INCREMENT SIZE
%NUMBER OF PARAMETER POINTS ON EACH AXIS [Delta_pb, T, R_0]
len       = zeros(1,3);
for k=1:3
    len(k) = uint8(round((U(k)-L(k))/increment(k))+1);
end

len       = [            21,    21,         21];

%Actual Parameter Space 
parameter_grid_delta = linspace(L(1,1), U(1,1),len(1)); 
parameter_grid_T     = linspace(L(1,2), U(1,2),len(2)); 
parameter_grid_R_0   = linspace(L(1,3), U(1,3),len(3));

%Chisquare Values over Parameter Space
chisquare_grid       = zeros(len(1), len(2), len(3));

disp(len);
for k = 1:len(1)%Vary Delta
    disp(k);
    for j = 1:len(2)%Vary T
        for l = 1:len(3)%Vary R_0
            chisquare_grid(k,j,l) = chisquare([parameter_grid_delta(k), parameter_grid_T(j), parameter_grid_R_0(l)]);
        end
    end
end


%Plotting Best Fit
%--------------------------------------------------------------------------
%Theoretical Fit
%Best Fit Paramaters
[minimum, min_index]              = min(chisquare_grid(:));
[delta_index, T_index, R_0_index] = ind2sub(size(chisquare_grid),min_index);
delta                             = parameter_grid_delta(delta_index);
T                                 = parameter_grid_T(T_index);
R_0                               = parameter_grid_R_0(R_0_index);

%Minimum Chisquare
[chisquare_min_val, itot_min] = SINcurr(delta, R_0, T, Input_V_j', Current_I_j', Total_Error_Current_I_j');

%Best Fit Parameter Estimate Uncertainties
Uncertainty_Index                 = zeros(2, 3); %min_index, max_index = index in parameter space where chi_square doubles for 3 parameters
Uncertainties                     = zeros(2, 3); %min,max = 2, numParameters = 3
Uncertainty_Avg                   = zeros(1,3);

%Delta Uncertainty
Uncertainty_Index(1,1)          = 1;
Uncertainty_Index(2,1)          = len(1);

for k= fliplr(1:delta_index)
    if chisquare_grid(k, T_index, R_0_index)/chisquare_min_val >= 2
        Uncertainty_Index(1,1)  = k;
        break
    end
end
for k = delta_index:len(1)
    if chisquare_grid(k, T_index, R_0_index)/chisquare_min_val >= 2
        Uncertainty_Index(2,1)  = k;
        break
    end
end

Uncertainties(1,1)  = parameter_grid_delta(Uncertainty_Index(1,1));
Uncertainties(2,1)  = parameter_grid_delta(Uncertainty_Index(2,1));
%T Uncertainty
Uncertainty_Index(1,2)          = 1;
Uncertainty_Index(2,2)          = len(2);
for k = fliplr(1:T_index)
    if chisquare_grid(delta_index, k, R_0_index)/chisquare_min_val >= 2
        Uncertainty_Index(1,2)  = k;
        break
    end
end
for k = T_index:len(2)
    if chisquare_grid(delta_index, k, R_0_index)/chisquare_min_val >= 2
        Uncertainty_Index(2,2)  = k;
        break
    end
end


Uncertainties(1,2)  = parameter_grid_T(Uncertainty_Index(1,2));
Uncertainties(2,2)  = parameter_grid_T(Uncertainty_Index(2,2));

%R_0 Uncertainty
Uncertainty_Index(1,3)          = 1;
Uncertainty_Index(2,3)          = len(3);
for k = fliplr(1:R_0_index)
    if chisquare_grid(delta_index, T_index, k)/chisquare_min_val >= 2
        Uncertainty_Index(1,3)  = k;
        break
    end
end
for k = R_0_index:len(3)
    if chisquare_grid(delta_index, T_index, k)/chisquare_min_val >= 2
        Uncertainty_Index(2,3)  = k;
        break
    end
end

Uncertainties(1,3)  = parameter_grid_R_0(Uncertainty_Index(1,3));
Uncertainties(2,3)  = parameter_grid_R_0(Uncertainty_Index(2,3));


%Plotting and Fitting
%--------------------------------------------------------------------------
%Linear Fit to find Approximate Resistance
trial_length = 255;%Number of Points in each Trial
i  = length(Input_V_j)/trial_length;%Number of Trials in each File
j  = 1;
colors = ['g','y','g','y'];
while j <= i
    x  = Input_V_j((j-1)*trial_length+1:j*trial_length,1);
    y  = Current_I_j((j-1)*trial_length+1:j*trial_length,1);
    dy = Total_Error_Current_I_j((j-1)*trial_length+1:j*trial_length,1);
    patch([x;flipud(x)],[y-dy;flipud(y+dy)],colors(j))
    hold on;
    j  = j + 1;
end
alpha(.01);

x  = Input_V_j;
y  = Current_I_j;
scatter(x,y, 'r.')
hold on;
scatter(Input_V_j, itot_min, 'g.')
hold on;
%axis([-inf inf -inf inf]);
axis([0 inf 0 inf]);
xlabel('Junction Voltage (V)');
ylabel('Junction Current (A)');
title(strcat(tunneling_type, ' I(V) Curve, Junction ', junction_type));


%Linear Fit Equation
Fit = polyfit(Input_V_j, Current_I_j,1);
equation = sprintf('', Fit(1), Fit(2));
%text(.001, 0, equation, 'FontSize', 10);
yL=get(gca,'YLim'); 
xL=get(gca,'XLim');   
%text((xL(1)+xL(2))/3,yL(2)*4/5,equation,...
text(0.02,0.0003,equation,...
      'HorizontalAlignment','left',...
      'VerticalAlignment','top',...
      'BackgroundColor',[1 1 1],...
      'FontSize',12);
plot(Input_V_j, Fit(1)*Input_V_j + Fit(2), 'r-.');



disp(chisquare_min_val);
disp(delta);
disp(T);
disp(R_0);

print(file_write, '-dpng');

function [chisquare_val] = chisquare(x)
    global Input_V_j;
    global Current_I_j;
    global Total_Error_Current_I_j;
    delta    = x(1,1);%meV
    T        = x(1,2);%meV
    R_0      = x(1,3);%Ohms
    Voltages = Input_V_j;
    Currents = Current_I_j;
    Error_Currents = Total_Error_Current_I_j;
    [chisquare_val, itot] = SINcurr(delta, R_0, T, Voltages', Currents', Error_Currents');
end

