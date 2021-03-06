clear all;
close all;
%Setup File Name
clear;
date               = '20190212'; %Date
date_written       = '20190218'; %Date Fitted
tunneling_type     =      'SIN'; %1: SIN, 2: NIN, 3: SIS
junction_type      =        '2'; %1,2, or 3
trial              =        '0'; %Just so we can cleanly store data, positive trials are forward sweeps, negative trials are backwards sweeps

%File Names
file_read          = strcat('measurementsAnalysis/CombinedErrorBars/', date, '_', junction_type, '_', tunneling_type, '_', 'Trial', string(trial), '.csv');
file_write         = strcat('FitFigures/', date_written, '_', junction_type, '_', tunneling_type, '_', 'Trial', string(trial), '.csv');

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
%Calculating R_0 Estimate
%Linear Fit Equation
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
    
disp(R_0_fit);

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


%Initialize
%[delta_pb (meV), delta_al (meV), T (K), R_0 (Ohms(] in mV/meV
L         = [.9,  .03,   1, R_0_fit-2];%LOWER BOUND
U         = [1.1,  .1,   2, R_0_fit+2];%UPPER BOUND
increment = [.05, .05, .05,        .2];%INCREMENT SIZE
%NUMBER OF PARAMETER POINTS ON EACH AXIS [Delta, T, R_0]
len       = zeros(1,4);
for k=1:4
    len(k) = uint8(round((U(k)-L(k))/increment(k))+1);
end

%Actual Parameter Space 
parameter_grid_delta_pb = linspace(L(1), U(1),len(1)); 
parameter_grid_delta_al = linspace(L(2), U(2),len(2)); 
parameter_grid_T        = linspace(L(3), U(3),len(3)); 
parameter_grid_R_0      = linspace(L(4), U(4),len(4));


%Chisquare Values over Parameter Space
chisquare_grid = zeros(len(1), len(2), len(3), len(4)); 

disp(len);

for k = 1:len(1)%Vary Delta_pb
    disp(k);
    for m = 1:len(2)%Vary Delta_al
        for j = 1:len(3)%Vary T
            for l = 1:len(4)%Vary R_0
                chisquare_grid(k, m, j, l) = chisquare([parameter_grid_delta_pb(k), parameter_grid_delta_al(m), parameter_grid_T(j), parameter_grid_R_0(l)]);
            end
        end
    end
end



%Theoretical Fit
%Fit Paramaterss
[minimum, min_index]              = min(chisquare_grid(:));
[delta_index, T_index, R_0_index] = ind2sub(size(chisquare_grid),min_index);
delta                             = parameter_grid_delta(delta_index);
T                                 = parameter_grid_T(T_index);
R_0                               = parameter_grid_R_0(R_0_index);


%Minimum Chisquare
[chisquare_min_val, itot_min] = SINcurr(delta, R_0, T, Input_V_j', Current_I_j', Total_Error_Current_I_j');

scatter(Input_V_j, Current_I_j, 'r.')
hold on;
scatter(Input_V_j, itot_min, 'g.')
disp(chisquare_min_val);
disp(delta);
disp(T);
disp(R_0);



function [chisquare_val] = chisquare(x)
    global Input_V_j;
    global Current_I_j;
    global Total_Error_Current_I_j;
    delta_pb = x(1,1);%meV
    delta_al = x(1,2);%meV
    T        = x(1,3);%K
    R_0      = x(1,4);%Ohms
    Voltages = Input_V_j;
    Currents = Current_I_j;
    Error_Currents = Total_Error_Current_I_j;
    [chisquare_val, itot] = SIScurr(delta_pb, delta_al, R_0, T, Voltages', Currents', Error_Currents');
end