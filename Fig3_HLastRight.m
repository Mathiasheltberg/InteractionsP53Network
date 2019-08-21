clear all; close all; clc
MM = 15;
Slopes = zeros(MM,2);
Amps = zeros(MM,4);
Amps2 = zeros(MM,4);
Pers = zeros(MM,3);
%  figure
%  set(gcf, 'Position', [744  830  240  191], 'Color', 'w');
%  set(gca, 'FontSize', 7);
for trial = 1:3
    trial
    for test = 1:MM
        Tmax = 600; ts = 0.1; dt = 0.001; N = 10; NN = round(Tmax/dt);
        Sav = zeros(round(Tmax/ts),3);
        P_I = ones(NN,1)*0.551;       M = ones(NN,1)*0.1;
        
        k1 = 0.1; TDel1 = round(2000);
        k2 = 1.0; k4 = 0.01;
        k8 = 0.1;  k7 = 0.15;
        
        k10 = k1; k20 = k2; k80 = k8; k70 = k7; k40 = k4;
        TDel10 = TDel1;
        
        RT = 0; click = 1; RT0 = 0; c = 1;
        X = 0;
        
        lambda1 = 3.;    lambda2 = 0.;    lambda3 = 0.15;
        rho1 = 1.0; rho2 = 1.0; kap1 = 1; kap2 = 1;
        if (trial == 1)
            rho1 = 1.0;
            kap2 = 0;
        elseif (trial == 2)
            rho2 = 0.01;
            kap1 = 0;
        elseif (trial == 3)
            rho2 = 0.1;
            kap1 = 0;
        end
        
        
        Timp = 300.0;
        k2 = k20*(1 + lambda1*X);
        k4 = k40*(1 + lambda2*X);
        k7 = k70*(1 - lambda3*X);
        k2;
        AR0 = 0; RTA = 0; ATR = 0;
        Peaks = []; CP = 0; Bot = []; CB = 0;
        SmPeaks = []; SmCP = 0; SmBot = []; SmCB = 0;
        TP = []; TB = []; SmTP = []; SmTB = [];
        
        while (RT < Tmax-dt)
            c = c+1;        
            RT  = RT + dt;                        
            
            if (RT > Timp)
                X = 0;
                if (RT < Timp + 10)
                    theta = 2*sqrt(test);
                else
                    theta = 0.1;
                end
                                
                k2 = k20 * (1 - kap1*ATR/(rho1 + ATR));
                ATR;
                k4 = k40 * (1 - kap2*ATR/(rho2 + ATR));
                ATR = ATR + dt*(theta - 1.0*(ATR));
            end
            %%%%%%%%%%%% Updates
            dPi = k1 - k2*M(c-1)*P_I(c-1)/(k4 + P_I(c-1));
            dM  = k7*(P_I(mod(c-TDel1,NN)+1)) - k8*M(c-1);
            
            P_I(c) = P_I(c-1) + dt*(dPi); M(c) = M(c-1) + dt*(dM);
            %%%%%%%%%%%%% Make sudden impace
            if (RT > Timp && RT0 < Timp);
                RTA = RT; click0 = click;
            end
            
            if (RT < 300 && RT > 200.0)
                if (P_I(c) < P_I(c-1) && P_I(c-2) < P_I(c-1))
                    SmCP = SmCP +1; SmPeaks(SmCP) =  P_I(c-1); SmTP(SmCP) = RT;
                end
                if (P_I(c) > P_I(c-1) && P_I(c-2) > P_I(c-1))
                    SmCB = SmCB +1; SmBot(SmCB) =  P_I(c-1); SmTB(SmCB) = RT;
                end
            end
            
            if (RT > 350 && RT < 480)
                if (P_I(c) < P_I(c-1) && P_I(c-2) < P_I(c-1))
                    CP = CP +1; Peaks(CP) =  P_I(c-1); TP(CP) = RT;
                end
                if (P_I(c) > P_I(c-1) && P_I(c-2) > P_I(c-1))
                    CB = CB +1; Bot(CB) =  P_I(c-1); TB(CB) = RT;
                end
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
%         figure(1)
%         plot(Sav(:,1),Sav(:,2)); hold on
%         plot(TP,Peaks,'o');

        AA = find(Sav(click0:click0+200,2)==max(Sav(click0:click0+200,2)));
        Amps2(test,trial) = mean((SmPeaks(end-3:end-1)-SmBot(end-3:end-1))/2);
        Amps(test,trial) = mean((Peaks(end-3:end-1)-Bot(end-3:end-1))/2)/Amps2(test,trial);
        Pers(test,trial) = mean(TP(2:end)-TP(1:end-1));
    end

%      set(gca, 'FontSize', 7);
%     plot(Amps(:,trial),'-'); hold on
     %   figure(2)
         set(gca, 'FontSize', 7);
         plot(Pers(:,trial)/18.3,'-'); hold on
         axis([0 15 0.0 2])
end


