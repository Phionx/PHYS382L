function chisquare_val = chisquare(x)
    %R_0 = ;
    %T = ;
    %Voltages =;
    %Currents = ;
    %Error_Currents = ;
    chisquare_val = SINcurr(x(1,1), x(1,2), R_0, T, Voltages, Currents, Error_Currents);
end