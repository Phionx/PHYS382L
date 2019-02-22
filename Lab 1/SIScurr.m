function [chisquare, itot] = SIScurr(delta_al, r_0, t, delta_pb, V_j, I_j, dI_j)
% energy and temperature measured in voltage units
%    1 meV = 1 mV
%    1 K   = 0.0861 mV

%delta_pb   = 2; %meV 

%Convert Temperature from K to meV
t  = t*.0861;


%Convert Voltages in V to mV
V_j = V_j*10^(3);
I_j = I_j*10^(3);
dI_j = dI_j*10^(3);


% Integration
e_max = 5; %cutoff for approximating integral
e_step = .01; %step size for approximating integral
e=-e_max:e_step:e_max;
[col,nptse]=size(e);

%Import voltage data
v = V_j;
[col,nptsv]=size(v);
%itot       = zeros(nptsv, 1);

% step over voltage
for k=1:nptsv
    vv=v(k);
    % do the integral over e
    for j=1:nptse
     if abs(e(j)+vv)>delta_pb
         % tame divergence in bcs density of states: .99 instead of 1
         ds_pb(j)=abs(e(j))/sqrt((e(j)+vv)^2-.99*delta_pb^2);
     else
         ds_pb(j)=0;
     end
     if abs(e(j))>delta_al
         % tame divergence in bcs density of states
         ds_al(j)=abs(e(j)+vv)/sqrt((e(j))^2-.99*delta_al^2);
     else
         ds_al(j)=0;
     end
     di(j)=ds_pb(j)*ds_al(j)*(1/(exp(e(j)/t)+1)-1/(exp((e(j)+vv)/t)+1));
    end
    % integrate with trapezoid rule include a parallel resistance r_par
    itot(k)=trapz(e,di)/r_0; %+ vv/r_par;
end

%Reduced Chi Square Calculation
%--------------------------------------------------------------------------
chisquare = 0;

for k =1:nptsv
    chisquare = chisquare + (I_j(k) - itot(k))^2/((dI_j(k))^2);
end

%nptsv = nptsv*(1000/100);%N*(f_filter)/(f_sample) accounting for freq filter correlation
numParameters = 3;
chisquare = (1/(nptsv-numParameters))*chisquare;

%Convert Current in mA to A
itot = itot*10^(-3);

end