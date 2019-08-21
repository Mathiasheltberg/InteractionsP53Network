clear all; close all; clc
Slopes = zeros(10,2);
Heights = zeros(10,2);
for trial = 1:2
    L1 = 1; L2 = 20;
    for test = L1:L2
        Tmax = 500; ts = 0.1; dt = 0.001; N = 10; NN = round(Tmax/dt);
        Sav = zeros(round(Tmax/ts),3);
        P_I = ones(NN,1)*0.551;       M = ones(NN,1)*0.1;
        
        k1 = 0.1; TDel1 = round(2000);
        k2 = 1.0; k4 = 0.01;
        k8 = 0.1;
        k7 = 0.15;
        k10 = k1; k20 = k2; k80 = k8; k70 = k7; k40 = k4;
        TDel10 = TDel1;
        
        RT = 0; click = 1; RT0 = 0; c = 1;
        X = 1;
        if (trial==1 || trial == 3 || trial == 5);  X = 0;  end
        
        lambda1 = 3.;    lambda2 = 0.;    lambda3 = 0.15;
        
        rho1 = 0.5;
        kap1 = 1; kap2 = 1;
        if (trial == 1 || trial == 2)
            kap2 = 0;
        elseif (trial == 3 || trial == 4)
            kap1 = 0;
        end
        Timp = 193.0;
        k2 = k20*(1 + lambda1*X);
        k4 = k40*(1 + lambda2*X);
        k7 = k70*(1 - lambda3*X);
        ATR = 0.0;
        AR0 = 0; RTA = 0; AR = 0;
        while (RT < Tmax-dt)
            c = c+1;        RT  = RT + dt;
            %k2 = k20 * (1 + lambda1*X) * (rho1/(rho1 + AR));
            if (RT > Timp)                
                if (RT < Timp + 10)
                    ATR = sqrt(4.0*test);
                    if (X < 0.5)
                        k2 = k20/(1+ATR);
                    elseif (X > 0.5)
                        k2 = k20*(1 + lambda1*X)/(1+ATR);
                    end
                else
                    ATR = ATR + dt*(0.1 - 1.0*(ATR));
                    k2 = k20*(1 + lambda1*X)/(1.0+ATR);
                end
            end
            %%%%%%%%%%%% Updates
            dPi = k1 - k2*M(c-1)*P_I(c-1)/(k4 + P_I(c-1));
            dM  = k7*(P_I(mod(c-TDel1,NN)+1)) - k8*M(c-1);
            
            P_I(c) = P_I(c-1) + dt*(dPi); M(c) = M(c-1) + dt*(dM);
            %%%%%%%%%%%%% Make sudden impace
            if (RT > Timp && RT0 < Timp);
                AR = 0.1;
                RTA = RT; click0 = click;
            end
            
            
            %%%%%%%%%%%% Save elements
            if (RT > ts*click);
                Sav(click,1) = RT;
                Sav(click,2) = P_I(c);
                Sav(click,4) = k2;
                Sav(click,3) = ATR;
                click=click+1; end
            RT0 = RT;
            
        end
        AA = find(Sav(click0:click0+400,2)==max(Sav(click0:click0+400,2)));
        plot(Sav(100:end-1,1),Sav(100:end-1,2)); hold on
        plot(Sav(click0+50,1),Sav(click0+50,2),'*')
        plot(Sav(click0,1),Sav(click0,2),'d')
        Slopes(test,trial) = (Sav(AA+click0,2) - Sav(click0,2))/(AA);
        Heights(test,trial) = (Sav(AA+click0,2));
        
        
    end
    figure; hold on
    plot(Sav(:,4),'d')
    plot(Sav(:,3),'*')
end
figure
%subplot(1,2,1)
plot(linspace(L1,L2,L2-L1+1)',Slopes(L1:L2,1),'b-*','Markersize',10,'LineWidth',2); hold on
plot(linspace(L1,L2,L2-L1+1)',Slopes(L1:L2,2),'r-*','Markersize',10,'LineWidth',2); goodplot
legend('Mdmx Supp','Control')
axis([0 20 Slopes(L1,2)/2 Slopes(L2,1)*2])
ylim([0 0.0085])
title('slope')
%goodplot
%subplot(1,2,2)
figure
plot(linspace(L1,L2,L2-L1+1)',Heights(L1:L2,1),'b-*','Markersize',10,'LineWidth',2); hold on
plot(linspace(L1,L2,L2-L1+1)',Heights(L1:L2,2),'r-*','Markersize',10,'LineWidth',2); goodplot
legend('Mdmx Supp','Control')
axis([0 20 Heights(L1,2)/2 Heights(L2,1)*2])
title('amplitude')
ylim([0 1])


