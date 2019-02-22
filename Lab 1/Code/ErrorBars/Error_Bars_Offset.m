%This matlab file is to combine the error bars in I_j and the error bars in
%V_j into combined error bars in I_j. Namely, we will propogate error in
%V_j using the approximate relation V_j/R_0 = I_j. Then we will add in 
%quadrature the original error in I_j with the error propagated from V_j 
%to get the total error in I_j. 
%Total_Error_I_j = \sqrt(Error_I_j^2 + \sqrt(


%Setup File Name
date_acquired      = '20190211'; %Date of Data Acquisition
date_written       = '20190212'; %Date of Data Analysis
tunneling_type     = 'SIS';%1: SIN, 2: NIN, 3: SIS
junction_type      = '3';%1,2, or 3
trial              = '-2';%Made so we can cleanly store data, positive trials are forward sweeps, negative trials are backwards sweep


file_reading              = strcat('measurementsAnalysis/', date_acquired, '_', junction_type, '_', tunneling_type, '_', 'Trial', string(trial), '.csv');
data                      = csvread(file_reading);
Current_I_j               = data(:, 1)';
Error_Current_I_j         = data(:, 2)';
Input_V_j                 = data(:, 3)';
Error_Input_V_j           = data(:,4)';

file_writing              = strcat('measurementsAnalysis/CorrectMeasurements/', date_written , '_', junction_type, '_', tunneling_type, '_', 'Trial', string(trial), '.csv');


%Initialization
%--------------------------------------------------------------------------
measurement_length = uint8(length(Current_I_j));%step number    
middle             = uint8(length(Current_I_j)/2);%middle
Offset             = Input_V_j(middle);
Error_Offset       = .001*Offset;
Input_V_j          = Input_V_j - Offset;

for k=1:measurement_length
    Error_Input_V_j(k) = sqrt((Error_Offset)^2 + (Error_Input_V_j(k))^2);
end

    
%Data Export
%--------------------------------------------------------------------------
export_data = [Current_I_j; Error_Current_I_j; Input_V_j; Error_Input_V_j];
export_data = export_data';

%Export to Right File
dlmwrite(file_writing, export_data, 'delimiter', ',', '-append');

