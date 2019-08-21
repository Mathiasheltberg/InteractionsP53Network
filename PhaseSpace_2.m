clear all; close all; clc

L = 10;
Val = zeros(100,100);
for test1 = 1:4
    test1
    for test2 = 1:1
        
    Tmax = 800; ts = 0.1; dt = 0.001; N = 10; NN = round(Tmax/dt);
    Sav = zeros(round(Tmax/ts),3);
    P_I = ones(NN,1)*0.1;       M = ones(NN,1)*0.15;
    d1 = 0; d2 = 0; d3 = 0; d4 = 0; d5 = 0; d6 = 0;
    k1 = 0.1; TDel1 = round(2200);
    k2 = 1.0; k4 = 0.01;
    k7 = 0.1; k8 = 0.1;
    k10 = k1; k20 = k2; k80 = k8; k70 = k7; k40 = k4;
    TDel10 = TDel1;

    RT = 0; click = 1; RT0 = 0; c = 1; X = 1.;
        
        lambda1 = 0; lambda2 = 0; lambda3 = 0; lambda4 = 0;
        lambda1 = test1/3;
        lambda2 = 0;
        k70 = 0.13;
        k20 = k20*(1 + lambda1*X);
        k40 = k40*(1 + lambda2*X);
        k70 = k70*(1 - lambda3*X);
        %k80 = k80*(1 - lambda4*X);
        TDel10 = TDel10*(1+lambda4*X);
        k1 = k10; k2 = k20; k4 = k40;
        k7 = k70; k8 = k80; TDel1 = TDel10;
        Tmax = 800;
        while (RT < Tmax-dt)
            c = c+1;        RT  = RT + dt;
            %%%%%%%%%%%% Updates
            
            dPi = k1 - k2*M(c-1)*P_I(c-1)/(k4 + P_I(c-1));
            dM  = k7*(P_I(mod(c-TDel1,NN)+1)) - k8*M(c-1);
            
            P_I(c) = P_I(c-1) + dt*(dPi);
            M(c)   = M(c-1) + dt*(dM);
            
            
            
            if (RT > 400 && RT0 < 400);
            k2 = k20/(1 + lambda1*X); k4 = k40/(1 + lambda2*X); 
            k7 = k70/(1 - lambda3*X); 
            TDel1 = TDel10/(1 + lambda4*X);
            click0 = click; 
            end
            
            %%%%%%%%%%%% Secure Against Runaways
            if (P_I(c) < 0);          P_I(c) = 0.000001;   end
            if (M(c) < 0);            M(c)   = 0.000001;   end
            RT0 = RT;
            %%%%%%%%%%%% Save elements
            if (RT > ts*click)
                SavDet(click,1) = RT; SavDet(click,2) = P_I(c); SavDet(click,3) = M(c); SavDet(click,4) = dPi; SavDet(click,5) = dM;
                click=click+1;
            end
        end
        AA = find(max(SavDet(click0:click0+1000,2)));
        Val(test1,test2) = SavDet(AA,2);
        if (test1 < 5)
            figure(1)
            plot(SavDet(click0-500:click0+1500,2),SavDet(click0-500:click0+1500,3),'LineWidth',2); hold on
            axis([0 0.3 0 0.3])
            pk = linspace(0.001,0.3,1000)';
            plot(pk,0.1*(1+0.01./pk),'--k')
            figure(2)
            plot(SavDet(click0-500:click0+1500,4),SavDet(click0-500:click0+1500,5),'LineWidth',2); hold on
            
        else
            plot(SavDet(click0-500:click0+1500,2),SavDet(click0-500:click0+1500,3)); hold on
        end
        
    end
end
xlabel('P_{53}'); ylabel('MDM2'); goodplot
