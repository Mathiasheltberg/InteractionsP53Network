clear all; close all; clc

for test = 1:3
    Tmax = 400; ts = 0.1; dt = 0.001; N = 10; NN = round(Tmax/dt);
    Sav = zeros(round(Tmax/ts),3);
    P_I = ones(NN,1)*0.551;       M = ones(NN,1)*0.1;
    d1 = 0; d2 = 0; d3 = 0; d4 = 0; d5 = 0; d6 = 0;
    k1 = 0.1; TDel1 = round(2000);
    k2 = 1.0; k4 = 0.01;
    k8 = 0.1;
    k7 = 0.15;
    k10 = k1; k20 = k2; k80 = k8; k70 = k7; k40 = k4;
    TDel10 = TDel1;
    
    RT = 0; click = 1; RT0 = 0; c = 1;
    
    X = 1;
    
    lambda1 = 3.;    lambda2 = 0.;    lambda3 = 0.15;
    
    Timp = 190.0;
    k2 = k20*(1 + lambda1*X);
    k4 = k40*(1 + lambda2*X);
    k7 = k70*(1 - lambda3*X);
    
    AR0 = 0; RTA = 0; AR = 0;
    
    while (RT < Tmax-dt)
        c = c+1;        RT  = RT + dt;
        
        if (RT > Timp)
            
            if (RT < Timp + 10)
                ATR = sqrt(4.0 * 8.0);
                k4 = k40*(1 + lambda2*X)*(0.1 - 0.01*test^2)/((0.1 - 0.01*test^2) + ATR);
            else
                ATR = ATR + dt*(0.1 - 1.0*(ATR));
                k4 = k40*(1 + lambda2*X)*(0.1 - 0.01*test^2)/((0.1 - 0.01*test^2) + ATR);
            end
        end
        
        %%%%%%%%%%%% Updates
        dPi = k1 - k2*M(c-1)*P_I(c-1)/(k4 + P_I(c-1));
        dM  = k7*(P_I(mod(c-TDel1,NN)+1)) - k8*M(c-1);
        
        P_I(c) = P_I(c-1) + dt*(dPi); M(c) = M(c-1) + dt*(dM);
        %%%%%%%%%%%%% Make sudden impace
        
        %%%%%%%%%%%% Secure Against Runaways
        if (P_I(c) < 0); P_I(c) = 0.000001; end;
        if (M(c) < 0);   M(c) = 0.000001; end;
        
        %%%%%%%%%%%% Save elements
        if (RT > ts*click);
            Sav(click,1) = RT;
            Sav(click,2) = P_I(c);
            Sav(click,4) = dPi;
            Sav(click,3) = M(c);
            click=click+1;
        end
        
        RT0 = RT;
    end
    plot(Sav(1:end-1,1),Sav(1:end-1,2),'LineWidth',3); hold on
    goodplot; axis([100 300 0 1.])
    
end


