%Setup File Name
date               = '20190212'; %Date of Data Acquisition
tunneling_type     = 'SIS';%1: SIN, 2: NIN, 3: SIS
junction_type      = '2';%1,2, or 3
trial              = '1';%Just so we can cleanly store data, positive trials are forward sweeps, negative trials are backwards sweeps
date_write         = '20190212'; %Writing Date
pic                = 'Fit0'; %convenient

file               = strcat('measurementsAnalysis/CombinedErrorBars/', date, '_', junction_type, '_', tunneling_type, '_', 'Trial', string(trial), '.csv');
file_write         = strcat('figures/', date_write, '_Junction', junction_type, '_', tunneling_type, '_', 'FIG', pic);
data = csvread(file);
Input_V_j = data(:, 1);
Current_I_j = data(:, 2);
Total_Error_Current_I_j = data(:, 3);




%Plotting and Fitting
%--------------------------------------------------------------------------
%Linear Fit to find Approximate Resistance
trial_length = 255;%Number of Points in each Trial
i  = length(Input_V_j)/trial_length;%Number of Trials in each File
j  = 1;
colors = ['b','y','b','y'];
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
%axis([-inf inf -inf inf]);
axis([0 inf 0 inf]);
xlabel('Junction Voltage (V)');
ylabel('Junction Current (A)');
title(strcat(tunneling_type, ' I(V) Curve, Junction ', junction_type));


%Linear Fit Equation
Fit = polyfit(Input_V_j, Current_I_j,1);
equation = sprintf('I = %.6f V + %.6f', Fit(1), Fit(2));
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
%print(file_write, '-dpng');
