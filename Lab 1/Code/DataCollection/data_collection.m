%Shantanu Jha & Chris West
%PHYS 382L: Superconducting Tunnel Junctions and the Josephson Effect
%SPRING 2019, Lab 1

%Signal Calibration
%--------------------------------------------------------------------------
AI            = analoginput('nidaq','dev1');
DIO           = digitalio('nidaq', 'dev1');
AO            = analogoutput('nidaq','dev1');
output        = addchannel(AO,1);
intput        = addchannel(AI,0);
step_duration = .01;%in seconds

%Set Sample Number per Step
set(AI, 'SampleRate',10000); %Rate for about 10000 samples/sec wanted
ActualRate    = get(AI,'SampleRate'); %Actual Rate Set
set(AI,'SamplesPerTrigger',step_duration*ActualRate); %Number of Samples/Step
set(AI,'TriggerType','Manual');

%Settings
%--------------------------------------------------------------------------
%file name
tunneling_type     = 'SIS';%1: SIN,2: NIN,3: SIS
junction_type      = '3';%1,2, or 3
trial              = '-1';%Just so we can cleanly store data
file               = strcat('Junction', junction_type, '_', tunneling_type, '_', 'Trial', string(trial), '.csv');


%data collection
direction          = -1;%+1 for forward sweep and -1 for backwards sweep
max_V_Output       = 5;%maximum output voltage
min_V_Output       = -5;%minimum output voltage
increment          = .02;%step size in our output voltage
delay              = .1; %How long should we extend the wait time for adequate data collection?

%Data Analysis
precision          = .0001;%uncertainty of output voltage
R_v                = 10000; %Variable Resistance, Ohms
Error_R_v          = R_v*.02; %Error in Variable Resistance
gain               = 1000;%gain from amplifier
offset             = -.00167;%volts
Error_offset       = .00001;%error
%Frequency_Filter   = 100Hz (Low Pass 6) 

%Initialization
%--------------------------------------------------------------------------
measurement_length = uint8(round((max_V_Output-min_V_Output)/increment)+1);%step number    


if direction == 1
    Output_V_c         = linspace(min_V_Output, max_V_Output, measurement_length);
elseif direction == -1
    Output_V_c         = linspace(max_V_Output, min_V_Output, measurement_length);
end

Error_Output_V_c   = zeros(1,measurement_length) + precision;

%Use V_c (output voltage) to get I_j (current across junction) = I_tot
Current_I_j        = Output_V_c/R_v;
Error_Current_I_j  = zeros(1,measurement_length);

%Propogate Error to Current across Junction
i = 1;
while i <= measurement_length
    Error_Current_I_j(i)  = sqrt((1/R_v)^2*(Error_Output_V_c(i))^2 + (Output_V_c(i)/(R_v^2))^2*(Error_R_v)^2);
    i = i + 1;
end




%Array to store data collected on input voltage (voltage V_j across junction)
Input_V_j          = zeros(1,measurement_length);
Error_Input_V_j    = zeros(1,measurement_length);


%Data Acquisition
%--------------------------------------------------------------------------
i = 1;
putdata(AO, 0);

    
figure(2);
while i <= measurement_length
    %SET Ouput Voltage
    putdata(AO, Output_V_c(i));
    
    %Open Signals
    start(AO);
    start(AI);
    
    
    %Trigger AI (Input) to start collecting data
    trigger(AI); %start data collection
    wait(AI, step_duration + delay); %Allow time for data collection
    data                = getdata(AI); %Store Data in variable
    
    %Storing average value and error in Junction Voltage Data Array
    Input_V_j(i)        = mean(data)/gain - offset; %rescaling for gain provided by amplifier, then subtracting offset
    if i == measurement_length || i == (measurement_length - 1)
        plot(data)
        hold on
    end
    
    Error_Input_V_j(i)  = sqrt((std(data)/gain)^2 + (Error_offset)^2);%Adding std deviation and error in our offset measurement in quadrature
    
    %Close Signals
    stop(AO);
    stop(AI);
    
    i = i + 1;
end

    
%Data Export
%--------------------------------------------------------------------------
export_data = [Current_I_j; Error_Current_I_j; Input_V_j; Error_Input_V_j];
export_data = export_data';

%Export to Right File
dlmwrite(file, export_data, 'delimiter', ',', '-append');
