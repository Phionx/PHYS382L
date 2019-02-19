%%%%%%%%%%%%%%%%%%%%
%   Copyright Jeremy Tanlimco, 2018.
%   Do not use without permission.
%   Ryan Lim (BK '19), Andrew Lingenfelter (MC '19), Jeremy Tanlimco (ES '19)
%   Physics 382: Advanced Physics Laboratory
%   Experiment 8: Superconducting Tunnel Junction and the Josephson Effect
%   Acquire, plot, and log I-V (current-voltage) data for an SIS junction.
%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all
openDAQ = daqfind;
for i = 1:length(openDAQ)
delete(openDAQ(i));
end

%%%%%%%%%%%%%%%%%%%%
%   Initialize DAQ input and output devices
%%%%%%%%%%%%%%%%%%%%
ao = analogoutput('Nidaq','Dev1');
out_chan = addchannel(ao,1);
ai = analoginput('Nidaq','Dev1');
in_chan = addchannel(ai,0);

%%%%%%%%%%%%%%%%%%%%
%   Initialize external variables
%   gain: Low-noise preamplifier gain
%   resistor: black box resistance translating voltage output to current
%%%%%%%%%%%%%%%%%%%%
gain = 500;
resistor = 1000;
junction = 1;

%%%%%%%%%%%%%%%%%%%%
%   Initialize various data acquisition parameters
%   voltage_step_size: voltage output step size
%   voltage_limit: upper and lower bounds of voltage output
%   time_step: time spent on each data point
%%%%%%%%%%%%%%%%%%%%
voltage_step_size = 0.00005;
voltage_limit = 0.005;
numSteps = (voltage_limit *2) / voltage_step_size + 1;
time_step = 0.01;
set(ai,'SampleRate',10000);
actualSampleRate = get(ai,'SampleRate');
set(ai,'SamplesPerTrigger',actualSampleRate * time_step);
set(ai,'TriggerType','Manual');

%%%%%%%%%%%%%%%%%%%%
%	Initialize voltage output and voltage input (measured data) arrays
%   Forward voltage: from -voltage_limit to voltage_limit
%%%%%%%%%%%%%%%%%%%%
dac_voltage_output_forward = linspace(-voltage_limit, voltage_limit,numSteps);
voltage_input_forward_average = zeros(1, numSteps);
voltage_input_forward_stdev = zeros(1, numSteps);

%%%%%%%%%%%%%%%%%%%%
%   Loop over voltage output and collect data into arrays
%   Take both average and standard deviation within the time step
%%%%%%%%%%%%%%%%%%%%
for i = 1:length(dac_voltage_output_forward)
    putdata(ao,dac_voltage_output_forward(i));
    start(ao);
    start(ai);
    trigger(ai);
    wait(ai,time_step + 0.25);
    [data,time] = getdata(ai);
    voltage_input_forward_average(i) = mean(data);
    voltage_input_forward_stdev(i) = std(data);
    stop(ai);
    stop(ao);
end

%%%%%%%%%%%%%%%%%%%%
%   Rescale voltage output to current and voltage inputs to correct for gain
%   Plot the data with error bars
%%%%%%%%%%%%%%%%%%%%
current_output_forward = dac_voltage_output_forward / resistor;
rescaled_voltage_input_forward_average = voltage_input_forward_average / gain;
rescaled_voltage_input_forward_stdev = voltage_input_forward_stdev / gain;
errorbar(current_output_forward,rescaled_voltage_input_forward_average,rescaled_voltage_input_forward_stdev,'LineStyle','none');
xlabel('Current to Junction (A)');
ylabel('Voltage across Junction (V)');
title('SIS I-V Curve Up');

%%%%%%%%%%%%%%%%%%%%
% Save the plot in an image file with a timestamp
%%%%%%%%%%%%%%%%%%%%
timeStamp = clock;
fix(timeStamp);
truncated_image_file_forward = sprintf('%d%d%d-%d%d-SIS-IV-Forward-Plot-%d.png',timeStamp(1),timeStamp(2),timeStamp(3),timeStamp(4),timeStamp(5), junction);
image_file_forward = fullfile('IV_Plots',truncated_image_file_forward);
saveas(gcf,image_file_forward,'png');


%%%%%%%%%%%%%%%%%%%%
% Save the current output and voltage outputs in a text file with a timestamp
%%%%%%%%%%%%%%%%%%%%
truncated_file_name_forward = sprintf('%d%d%d-%d%d-SIS-IV-Forward-%d.txt',timeStamp(1),timeStamp(2),timeStamp(3),timeStamp(4),timeStamp(5), junction);
file_name_forward = fullfile('IV_Data',truncated_file_name_forward);
file_forward = fopen(file_name_forward,'w');
fprintf(file_forward,'%8s\t%8s\t%8s\t%8s\r\n','Gain','Resistance','Time Step','Junction');
fprintf(file_forward,'%8.6f\t%8.6f\t%8.6f\t%8.6f\r\n',[gain,resistor,time_step,junction]);
fprintf(file_forward,'%8s\t%8s\t%8s\r\n','I Out','Avg V Across', 'StdDev V Across');
fprintf(file_forward,'%8.6f\t%8.6f\t%8.6f\r\n',[current_output_forward; rescaled_voltage_input_forward_average; rescaled_voltage_input_forward_stdev]);
fclose(file_forward);

%%%%%%%%%%%%%%%%%%%% 
%	Initialize voltage output and voltage input (measured data) arrays
%   Backward voltage: from -voltage_limit to voltage_limit
%%%%%%%%%%%%%%%%%%%%
dac_voltage_output_backward = linspace(voltage_limit, -voltage_limit,numSteps);
voltage_input_backward_average = zeros(1, numSteps);
voltage_input_backward_stdev = zeros(1, numSteps);

%%%%%%%%%%%%%%%%%%%%
%   Loop over voltage output and collect data into arrays
%   Take both average and standard deviation within the time step
%%%%%%%%%%%%%%%%%%%%
for i = 1:length(dac_voltage_output_backward)
    putdata(ao,dac_voltage_output_backward(i));
    start(ao);
    start(ai);
    trigger(ai);
    wait(ai,time_step + 0.25);
    [data,time] = getdata(ai);
    voltage_input_backward_average(i) = mean(data);
    voltage_input_backward_stdev(i) = std(data);
    stop(ai);
    stop(ao);
end

%%%%%%%%%%%%%%%%%%%%
%   Rescale voltage output to current and voltage inputs to correct for gain
%   Plot the data with error bars
%%%%%%%%%%%%%%%%%%%%
current_output_backward = dac_voltage_output_backward / resistor;
rescaled_voltage_input_backward_average = voltage_input_backward_average / gain;
rescaled_voltage_input_backward_stdev = voltage_input_backward_stdev / gain;
errorbar(current_output_backward,rescaled_voltage_input_backward_average,rescaled_voltage_input_backward_stdev,'LineStyle','none');
xlabel('Current to Junction (A)');
ylabel('Voltage across Junction (V)');
title('SIS I-V Curve Down');

%%%%%%%%%%%%%%%%%%%%
% Save the plot in an image file with a timestamp
%%%%%%%%%%%%%%%%%%%%
timeStamp = clock;
fix(timeStamp);
truncated_image_file_backward = sprintf('%d%d%d-%d%d-SIS-IV-Backward-Plot-%d.png',timeStamp(1),timeStamp(2),timeStamp(3),timeStamp(4),timeStamp(5),junction);
image_file_backward = fullfile('IV_Plots',truncated_image_file_backward);
saveas(gcf,image_file_backward,'png');

%%%%%%%%%%%%%%%%%%%%
% Save the current output and voltage outputs in a text file with a timestamp
%%%%%%%%%%%%%%%%%%%%
truncated_file_name_backward = sprintf('%d%d%d-%d%d-SIS-IV-Backward-%d.txt',timeStamp(1),timeStamp(2),timeStamp(3),timeStamp(4),timeStamp(5),junction);
file_name_backward = fullfile('IV_Data',truncated_file_name_backward);
file_backward = fopen(file_name_backward,'w');
fprintf(file_backward,'%8s\t%8s\t%8s\t%8s\r\n','Gain','Resistance','Time Step','Junction');
fprintf(file_backward,'%8.6f\t%8.6f\t%8.6f\t%8.6f\r\n',[gain,resistor,time_step,junction]);
fprintf(file_backward,'%8s\t%8s\t%8s\r\n','I Out','Avg V Across', 'StdDev V Across');
fprintf(file_backward,'%8.6f\t%8.6f\t%8.6f\r\n',[current_output_backward; rescaled_voltage_input_backward_average; rescaled_voltage_input_backward_stdev]);
fclose(file_backward);
