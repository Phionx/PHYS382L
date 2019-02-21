%This matlab file is to get rid of our Current_I_j error bars because Prof
%Lamoreaux said they were unnecessary


%Setup File Name
date_acquired      = '20190212'; %Date of Data Acquisition
date_written       = '20190219'; %Date of Data Analysis
tunneling_type     = 'SIS';%1: SIN, 2: NIN, 3: SIS
junction_type      = '3';%1,2, or 3
trial              = '2';%Made so we can cleanly store data, positive trials are forward sweeps, negative trials are backwards sweeps


file_reading              = strcat('measurementsAnalysis/CombinedTrials/', date_acquired, '_', junction_type, '_', tunneling_type, '_', 'Trial', string(trial), '.csv');
data                      = csvread(file_reading);
Current_I_j               = data(:, 1)';
Input_V_j                 = data(:, 3)';
Error_Input_V_j           = data(:,4)';

file_writing              = strcat('measurementsAnalysis/LamoreauxErrorBars/', date_written , '_', junction_type, '_', tunneling_type, '_', 'Trial', string(trial), '.csv');


%Get Rid of useless Current Error Bars
%--------------------------------------------------------------------------
Error_Current_I_j         = zeros(1,length(Current_I_j));

%Data Export
%--------------------------------------------------------------------------
export_data = [Current_I_j; Error_Current_I_j; Input_V_j; Error_Input_V_j];
export_data = export_data';

%Export to Right File
dlmwrite(file_writing, export_data, 'delimiter', ',', '-append');

