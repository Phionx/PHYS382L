function [chisquare, itot] = NINcurr(r_0, V_j, I_j, dI_j)
% energy and temperature measured in voltage units
%    1 meV = 1 mV
%    1 K   = 0.0861 mV

%Convert Voltages in V to mV
V_j  = V_j*10^(3);
I_j  = I_j*10^(3);
dI_j = dI_j*10^(3);

itot = V_j*(1/r_0);

nptsv = length(V_j);
%Reduced Chi Square Calculation
%--------------------------------------------------------------------------
chisquare = 0;

for k =1:nptsv
    chisquare = chisquare + (I_j(k) - itot(k))^2/((dI_j(k))^2);
end

%nptsv = nptsv*(1000/100);%N*(f_filter)/(f_sample) accounting for freq filter correlation
numParameters = 1;
chisquare = (1/(nptsv-numParameters))*chisquare;

%Convert Current in mA to A
itot = itot*10^(-3);

end