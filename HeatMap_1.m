clear all; close all; clc

L = 40;
Val = zeros(L,L);
for test1 = 1:L
    test1
    t = linspace(1,20,40);
    for test2 = 1:3
        
        Tmax = 400; ts = 0.1; dt = 0.001; N = 10; NN = round(Tmax/dt);
        Sav = zeros(round(Tmax/ts),3);
        P_I = ones(NN,1)*0.1;       M = ones(NN,1)*0.15;
        d1 = 0; d2 = 0; d3 = 0; d4 = 0; d5 = 0; d6 = 0;
        k1 = 0.1; TDel1 = round(2200);
        k2 = 1.0; k4 = 0.01;
        k7 = 0.1; k8 = 0.1;
        k10 = k1; k20 = k2; k80 = k8; k70 = k7; k40 = k4;
        TDel10 = TDel1;
        
        RT = 0; click = 1; RT0 = 0; c = 1;
        X = 1.;
        
        lambda1 = 0; lambda2 = 0; lambda3 = 0; lambda4 = 0;
        lambda1 = test1/2;
        lambda3 = test2/4.0;
        k70 = 0.13;
        k20 = k20*(1 + lambda1*X);
        k40 = k40*(1 + lambda2*X);
        k70 = k70*(1 - lambda3*X);
        
        TDel10 = TDel10*(1+lambda4*X);
        k1 = k10; k2 = k20; k4 = k40;
        k7 = k70; k8 = k80; TDel1 = TDel10;
        
        while (RT < Tmax-dt)
            c = c+1;        RT  = RT + dt;
            %%%%%%%%%%%% Updates
            
            dPi = k1 - k2*M(c-1)*P_I(c-1)/(k4 + P_I(c-1));
            dM  = k7*(P_I(mod(c-TDel1,NN)+1)) - k8*M(c-1);
            
            P_I(c) = P_I(c-1) + dt*(dPi);
            M(c)   = M(c-1) + dt*(dM);
            
            
            
            if (RT > 150 && RT0 < 150);
                k2 = k20/(1 + lambda1*X); k4 = k40/(1 + lambda2*X);
                k7 = k70/(1 - lambda3*X);
                TDel1 = TDel10/(1 + lambda4*X);
                click0 = click; end
            
            %%%%%%%%%%%% Secure Against Runaways
            if (P_I(c) < 0);          P_I(c) = 0.000001;   end
            if (M(c) < 0);            M(c)   = 0.000001;   end
            RT0 = RT;
            %%%%%%%%%%%% Save elements
            if (RT > ts*click)
                SavDet(click,1) = RT; SavDet(click,2) = P_I(c); SavDet(click,3) = M(c);
                click=click+1;
            end
        end
        AA = find(SavDet(click0:click0+200,2)==max(SavDet(click0:click0+200,2)));
        Val(test1,test2) = SavDet(AA+click0,2);
        
        plot(SavDet(:,1),SavDet(:,2)); hold on
        plot(SavDet(AA+click0,1),SavDet(AA+click0,2),'*')
        
    end
    
end
xlabel('P_{53}'); ylabel('MDM2'); goodplot
figure

plot(t,Val(:,1),'*'); hold on; goodplot
plot(t,Val(:,2),'*'); hold on; goodplot
plot(t,Val(:,3),'*'); hold on; goodplot
plot(t,Val(:,4),'*'); hold on; goodplot
legend(texlabel('lambda_3 = +'),texlabel('lambda_3 = ++'),texlabel('lambda_3 = +++'))
figure

semilogy(t,Val(:,1),'*'); hold on; goodplot
semilogy(t,Val(:,2),'*'); hold on; goodplot
semilogy(t,Val(:,3),'*'); hold on; goodplot
semilogy(t,Val(:,4),'*'); hold on; goodplot
legend(texlabel('lambda_3 = +'),texlabel('lambda_3 = ++'),texlabel('lambda_3 = +++'))
% figure
% imagesc(Val); colormap hot; goodplot

