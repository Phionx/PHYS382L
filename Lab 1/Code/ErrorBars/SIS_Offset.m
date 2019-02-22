%This matlab file is to correct the asymmetry in our SIS Data


%Setup File Name
date_acquired      = '20190220'; %Date of Data Acquisition
date_written       = '20190221'; %Date of Data Analysis
tunneling_type     = 'SIS';%1: SIN, 2: NIN, 3: SIS
junction_type      = '3';%1,2, or 3
trial              = '0';%Made so we can cleanly store data, positive trials are forward sweeps, negative trials are backwards sweeps

file_reading              = strcat('Data/measurementsAnalysis/FinalData/', date_acquired, '_', junction_type, '_', tunneling_type, '_', 'Trial', string(trial), '.csv');
data                      = csvread(file_reading);
Input_V_j                 = data(:, 1)';
Current_I_j               = data(:, 2)';
Error_Current_I_j         = data(:, 3)';

file_writing              = strcat('Data/measurementsAnalysis/FinalDataSISOffset/', date_written , '_', junction_type, '_', tunneling_type, '_', 'Trial', string(trial), '.csv');


Input_V_j                 = Input_V_j + .0805*10^(-3);

    
%Data Export
%--------------------------------------------------------------------------
export_data = [Input_V_j; Current_I_j; Error_Current_I_j];
export_data = export_data';

%Export to Right File
dlmwrite(file_writing, export_data, 'delimiter', ',');

