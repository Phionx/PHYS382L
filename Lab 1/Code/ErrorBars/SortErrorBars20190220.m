%This matlab file is to sort our data points in increasing order of their
%Input_V_j values (x values) 


%Setup File Name
date_acquired      = '20190212'; %Date of Data Acquisition
date_written       = '20190221'; %Date of Data Analysis
tunneling_type     = 'SIS';%1: SIN, 2: NIN, 3: SIS
junction_type      = '3';%1,2, or 3
trial              = '0';%Made so we can cleanly store data, positive trials are forward sweeps, negative trials are backwards sweeps


file_reading              = strcat('Data/measurementsAnalysis/WrongAnalysis/CombinedErrorBars/', date_acquired, '_', junction_type, '_', tunneling_type, '_', 'Trial', string(trial), '.csv');
disp(file_reading);
data                      = csvread(file_reading);
data                      = sortrows(data); %sorts by first row values, this is to get our Input_V_j in order!
file_writing              = strcat('Data/measurementsAnalysis/SortedWrongAnalysis/', date_written , '_', junction_type, '_', tunneling_type, '_', 'Trial', string(trial), '.csv');

%Data Export
%--------------------------------------------------------------------------
export_data = data;

%Export to Right File
dlmwrite(file_writing, export_data, 'delimiter', ',', '-append');

