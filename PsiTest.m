clear all; close all; clc

DisSav = zeros(100,2);
for test = 1:3
    
    Tmax = 800; ts = 0.1; dt = 0.001; N = 10; NN = round(Tmax/dt);
    Sav = zeros(round(Tmax/ts),3);
    P_I = ones(NN,1)*0.1*test;       M = ones(NN,1)*0.15;
    d1 = 0; d2 = 0; d3 = 0; d4 = 0; d5 = 0; d6 = 0;
    k1 = 0.1; TDel1 = round(2000);
    k2 = 1.; k4 = 0.01;
    k7 = 0.15; k8 = 0.1;
    k10 = k1; k20 = k2; k80 = k8; k70 = k7; k40 = k4;
    TDel10 = TDel1;
    
    RT = 0; click = 1; RT0 = 0; c = 1;
    
    X = 1.;
    
    lambda1 = 0; lambda2 = 0; lambda3 = 0; lambda4 = 0;
    %%%%%%%% THIS IS THE LAMBDA YOU CHANGE
    lambda1 = 3.0*test;
    %   lambda2 = test*2.0;
    % lambda3 = test*0.25;
    %  lambda4 = test*1.0;
%     P_I(1) = 0.1*test;
    k20 = k20*(1 + lambda1*X);
    k40 = k40*(1 + lambda2*X);
    k70 = k70*(1 - lambda3*X);
    TDel10 = TDel10*(1+lambda4*X);
    
    k1 = k10; k2 = k20; k4 = k40;
    k7 = k70; k8 = k80; TDel1 = TDel10;
    click0 = 0
    while (RT < Tmax-dt)
        c = c+1;        RT  = RT + dt;
        %%%%%%%%%%%% Updates
        
        dPi = k1 - k2*M(c-1)*P_I(c-1)/(k4 + P_I(c-1));
        dM  = k7*(P_I(mod(c-TDel1,NN)+1)) - k8*M(c-1);
        
        P_I(c) = P_I(c-1) + dt*(dPi);
        M(c) = M(c-1) + dt*(dM);
        
        
        %%%%%%%%%%%%% Make sudden impace
        if (RT > 200 && RT0 < 200);
            k2 = k20/(1 + lambda1*X); k4 = k40/(1 + lambda2*X);
            k7 = k70/(1 - lambda3*X);
            % k8 = k80/(1 - lambda4*X);
            TDel1 = TDel10/(1 + lambda4*X);
            click0 = click; end
        
        %%%%%%%%%%%% Secure Against Runaways
        if (P_I(c) < 0); P_I(c) = 0.000001; end;
        if (M(c) < 0);   M(c) = 0.000001; end;
        
        %%%%%%%%%%%% Save elements
        if (RT > ts*click);
            Sav(click,1) = RT;
            Sav(click,2) = P_I(c);
            Sav(click,3) = M(c);
            click=click+1; end
        
        RT0 = RT;
        
    end
    
    plot(Sav(1:click0+500,1),Sav(1:click0+500,2),'LineWidth',2); hold on
    if (test == 3)
        plot(Sav(click0+600:end-1,1),Sav(click0+600:end-1,2),'k','LineWidth',2); hold on;
    end
    
end
axis([0 300 0 0.5])
xlabel('Time [AU]'); ylabel('P53 [AU]'); goodplot;
%legend(texlabel('lambda_3 = +'), texlabel('lambda_3 = ++'),texlabel('lambda_3 = +++'));
