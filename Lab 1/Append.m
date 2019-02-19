%This matlab file is to combine data sets from different trials into one
%file


%Setup File Name
date_source        = '20190212'; %Date of Data Acquisition of data to append
tunneling_type     = 'SIS';%1: SIN, 2: NIN, 3: SIS
junction_type      = '3';%1,2, or 3
trial              = '2';%Made so we can cleanly store data, positive trials are forward sweeps, negative trials are backwards sweeps

date_dest               = '20190212'; %Date of Data Acquisition of file in which data will be appended
tunneling_type_dest     = 'SIS';%1: SIN, 2: NIN, 3: SIS
junction_type_dest      = '3';%1,2, or 3
trial_dest              = '0';%Made so we can cleanly store data, positive trials are forward sweeps, negative trials are backwards sweeps


file_reading              = strcat('measurementsAnalysis/CombinedTrials/', date_source, '_', junction_type, '_', tunneling_type, '_', 'Trial', string(trial), '.csv');
data                      = csvread(file_reading);
Current_I_j               = data(:, 1)';
Error_Current_I_j         = data(:, 2)';
Input_V_j                 = data(:, 3)';
Error_Input_V_j           = data(:, 4)';

file_writing              = strcat('measurementsAnalysis/CombinedTrials/', date_dest , '_', junction_type_dest, '_', tunneling_type_dest, '_', 'Trial', string(trial_dest), '.csv');

%Data Write
%--------------------------------------------------------------------------
export_data = [Current_I_j; Error_Current_I_j; Input_V_j; Error_Input_V_j];
export_data = export_data';

%Export to Right File
dlmwrite(file_writing, export_data, 'delimiter', ',', '-append');

