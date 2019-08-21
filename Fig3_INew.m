clear all; close all; clc
Slopes = zeros(10,2);
Heights = zeros(10,2);
for trial = 1:2
    for test = 8:16
        Tmax = 400; ts = 0.1; dt = 0.001; N = 10; NN = round(Tmax/dt);
        Sav = zeros(round(Tmax/ts),3);
        P_I = ones(NN,1)*0.551;       M = ones(NN,1)*0.1;
        
        k1 = 0.1; TDel1 = round(2000);
        k2 = 1.0; k4 = 0.01;
        k8 = 0.1;
        k7 = 0.13;
        k10 = k1; k20 = k2; k80 = k8; k70 = k7; k40 = k4;
        TDel10 = TDel1;
        
        RT = 0; click = 1; RT0 = 0; c = 1;
        X = 1;
        if (trial==1 || trial == 3 || trial == 5);  X = 0;  end
        
        lambda1 = 3.;    lambda2 = 0.;    lambda3 = 0.15;
        
        rho1 = 0.1; 
        kap1 = 1; kap2 = 1;
        if (trial == 1 || trial == 2)
            kap2 = 0;
        elseif (trial == 3 || trial == 4)
            kap1 = 0;
        end
        Timp = 180.0;
        k2 = k20*(1 + lambda1*X);
        k4 = k40*(1 + lambda2*X);
        k7 = k70*(1 - lambda3*X);
        
        AR0 = 0; RTA = 0; AR = 0;
        rho1 = rho1*(1.+2*X);
        while (RT < Tmax-dt)
            c = c+1;        RT  = RT + dt;
            k2 = k20 * (1 + lambda1*X) * (rho1/(rho1 + AR));
            
            %%%%%%%%%%%% Updates
            dPi = k1 - k2*M(c-1)*P_I(c-1)/(k4 + P_I(c-1));
            dM  = k7*(P_I(mod(c-TDel1,NN)+1)) - k8*M(c-1);
            
            P_I(c) = P_I(c-1) + dt*(dPi); M(c) = M(c-1) + dt*(dM);
            %%%%%%%%%%%%% Make sudden impace
            if (RT > Timp && RT0 < Timp);
                
                if (test == 0)
                    X = 0; k2 = k20*(1 + lambda1*X); k4 = k40*(1 + lambda2*X); k7 = k70*(1 - lambda3*X);
                else
                    AR = 0.1*test
                    RTA = RT; click0 = click;
                end
            end
            if (RT > Timp );
               if (RT < Timp + test)
                    AR = 0.1*test;
                else
                    AR = AR + dt*(0.01 - 0.1*AR);
               end
            end
            %%%%%%%%%%%% Secure Against Runaways
            if (P_I(c) < 0); P_I(c) = 0.000001; end;
            if (M(c) < 0);   M(c) = 0.000001; end;
            
            %%%%%%%%%%%% Save elements
            if (RT > ts*click);
                Sav(click,1) = RT;
                Sav(click,2) = P_I(c);
                Sav(click,4) = dPi;
                Sav(click,3) = M(c);
                click=click+1; end            
            RT0 = RT;
            
        end
        AA = find(Sav(click0:click0+200,2)==max(Sav(click0:click0+200,2)));
        plot(Sav(100:end-1,1),Sav(100:end-1,2)); hold on
        plot(Sav(click0+AA,1),Sav(click0+AA,2),'*')
        plot(Sav(click0+40,1),Sav(click0+40,2),'d')
        Slopes(test,trial) = mean(Sav(click0:click0+20,4));
        Heights(test,trial) = (Sav(AA+click0,2));
    end
end
figure
plot(linspace(8,16,9)',Slopes(8:16,1),'b-*','Markersize',10,'LineWidth',2); hold on
plot(linspace(8,16,9)',Slopes(8:16,2),'r-*','Markersize',10,'LineWidth',2); goodplot
legend('Mdmx Supp','Control')
axis([6 18 0 0.1])
figure
plot(linspace(8,16,9)',Heights(8:16,1),'b-*','Markersize',10,'LineWidth',2); hold on
plot(linspace(8,16,9)',Heights(8:16,2),'r-*','Markersize',10,'LineWidth',2); goodplot
legend('Mdmx Supp','Control')
axis([6 18 0 1.0])


