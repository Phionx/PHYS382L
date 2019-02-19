%This matlab file is to combine the error bars in I_j and the error bars in
%V_j into combined error bars in I_j. Namely, we will propogate error in
%V_j using the approximate relation V_j/R_0 = I_j. Then we will add in 
%quadrature the original error in I_j with the error propagated from V_j 
%to get the total error in I_j. 
%Total_Error_I_j = \sqrt(Error_I_j^2 + \sqrt(


%Setup File Name
date_acquired      = '20190205'; %Date of Data Acquisition
date_written       = '20190211'; %Date of Data Analysis
tunneling_type     = 'SISC';%1: SIN, 2: NIN, 3: SIS
junction_type      = '3';%1,2, or 3
trial              = '-1';%Made so we can cleanly store data, positive trials are forward sweeps, negative trials are backwards sweeps

%data collection
max_V_Output       = 5;%maximum output voltage
min_V_Output       = -5;%minimum output voltage
increment          = .02;%step size in our output voltage
direction          = -1;%direction of sweep

%Data Analysis
%CURRENT (I_j)
precision_output   = .005;%uncertainty of driving voltage (used to get current I_j)
R_v                = 100000; %Variable Resistance, Ohms
Error_R_v          = R_v*.02; %Error in Variable Resistance

%Voltage (V_j)
precision_input    = .005;%uncertainty of measured voltage (V_j)
gain               = 500;%gain from amplifier
offset             = -.00167;%volts
Error_offset       = .00001;%error


file_reading              = strcat('measurements/', date_acquired, '_', junction_type, '_', tunneling_type, '_', 'Trial', string(trial), '.csv');
data                      = csvread(file_reading);
Current_I_j               = data(:, 1)';
Input_V_j                 = data(:, 3)';
Error_Input_V_j           = data(:,4)';

file_writing              = strcat('measurements/CorrectErrorBars/', date_written , '_', junction_type, '_', tunneling_type, '_', 'Trial', string(trial), '.csv');


%Initialization
%--------------------------------------------------------------------------
measurement_length = uint8(round((max_V_Output-min_V_Output)/increment)+1);%step number    



%CURRENT ERROR BARS -------------------------------------------------------
if direction == 1
    Output_V_c         = linspace(min_V_Output, max_V_Output, measurement_length);
elseif direction == -1
    Output_V_c         = linspace(max_V_Output, min_V_Output, measurement_length);
end

Error_Output_V_c   = zeros(1,measurement_length) + precision_output;

%Use V_c (output voltage) to get I_j (current across junction) = I_tot
Error_Current_I_j  = zeros(1,measurement_length);

%Propogate Error to Current across Junction
i = 1;
while i <= measurement_length
    Error_Current_I_j(i)  = sqrt((1/R_v)^2*(Error_Output_V_c(i))^2 + (Output_V_c(i)/(R_v^2))^2*(Error_R_v)^2);
    i = i + 1;
end

    
%Data Export
%--------------------------------------------------------------------------
export_data = [Current_I_j; Error_Current_I_j; Input_V_j; Error_Input_V_j];
export_data = export_data';

%Export to Right File
dlmwrite(file_writing, export_data, 'delimiter', ',', '-append');

