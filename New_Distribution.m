clear all; close all; clc
MM = 1000;
DisSav = zeros(MM,4);
Tmax = 200; ts = 1; dt = 0.001; N = 10;
Pe = zeros(MM,1); PIAL = zeros(round(Tmax/ts),MM);
enois = 0.06; inois = 0.06;
startP =0.1320; startM = 0.1192;
for testTT = 1:3
    c = 1;
    figure; hold on
    for testN = 1:MM
        testN
        NN = round(Tmax/dt);
        Sav = zeros(round(Tmax/ts),3);
        P_I = ones(NN,1)*startP;       M = ones(NN,1)*startM;
        d1 = 0; d2 = 0; d3 = 0; d4 = 0; d5 = 0; d6 = 0;
        k1 = 0.1; TDel1 = round(2000);
        k2 = 1.0; k4 = 0.01;
        k7 = 0.13; k8 = 0.1;
        k10 = k1; k20 = k2; k80 = k8; k70 = k7; k40 = k4;
        TDel10 = TDel1;
        
        RT = 0; click = 1; RT0 = 0; c = 1;
        
        X = 1.;
        
        lambda1 = 3.; lambda2 = 0.; lambda3 = 0.15;
        if (testTT == 2)
            X = 0;
        elseif (testTT == 3)
            lambda3 = 0.0;
        elseif (testTT == 4)
            TDel10 = round(2000); X = 0.0;
        end
        k20 = k20*(1 + lambda1*X);
        k40 = k40*(1 + lambda2*X);
        k70 = k70*(1 - lambda3*X);
        if (testN < 1000)
            d1 = enois*randn; k10 = k10*(1 + d1);
            d2 = enois*randn; k20 = k20*(1 + d2);
            d3 = enois*randn; k40 = k40*(1 + d3);
            d4 = enois*randn; k70 = k70*(1 + d4);
            d5 = enois*randn; k80 = k80*(1 + d5);
            d6 = enois*randn; TDel1 = round(TDel1*(1 + d6));
        end
        k1 = k10; k2 = k20; k4 = k40; k7 = k70; k8 = k80; TDel1 = TDel10;
        
        
        
        while (RT < Tmax-dt)
            c = c+1;        RT  = RT + dt;
            %%%%%%%%%%%% Updates
            
            dPi = k1 - k2*M(c-1)*P_I(c-1)/(k4 + P_I(c-1));
            dM  = k7*(P_I(mod(c-TDel1,NN)+1)) - k8*M(c-1);
            
            P_I(c) = P_I(c-1) + dt*(dPi) + sqrt(dt)*(sqrt(k1) - sqrt(k2*M(c-1)*P_I(c-1)/(k4 + P_I(c-1))))*inois*randn;
            M(c) = M(c-1) + dt*(dM) + sqrt(dt)*(sqrt( k7*(P_I(mod(c-TDel1,NN)+1))) - sqrt(k8*M(c-1)))*inois*randn;
%             
%             P_I(c) = P_I(c-1) + dt*(dPi) + inois*randn;
%             M(c) = M(c-1) + dt*(dM) + inois*randn;
            
            %%%%%%%%%%%% Secure Against Runaways
            if (P_I(c) < 0); P_I(c) = 0.000001; end;
            if (M(c) < 0);   M(c) = 0.000001; end;
            
            if (RT > ts*click);
                Sav(click,1) = RT;
                Sav(click,2) = P_I(c);
                Sav(click,3) = M(c);
                click=click+1; end
            
            
        end        
        DisSav(testN,testTT) = M(c)/P_I(c);
        startP = P_I(c);
        startM = M(c);
        plot(Sav(1:end-1,1),Sav(1:end-1,2)); axis([0 100 0 0.4]); goodplot

    end
   
end

le = 0.0:0.02:2.5;
Dny = zeros(length(le)-1,3);
Dny(1,1) = length(find(DisSav(:,1) < le(1) ) );
Dny(1,2) = length(find(DisSav(:,2) < le(1) ) );
Dny(1,3) = length(find(DisSav(:,3) < le(1) ) );
%Dny(1,4) = length(find(DisSav(:,4) < le(1) ) );
for i = 2:length(le)-1;
    Dny(i,1) = length(find(DisSav(:,1) < le(i) & DisSav(:,1) >= le(i-1) ));
    Dny(i,2) = length(find(DisSav(:,2) < le(i) & DisSav(:,2) >= le(i-1) ));
    Dny(i,3) = length(find(DisSav(:,3) < le(i) & DisSav(:,3) >= le(i-1) ));
   % Dny(i,4) = length(find(DisSav(:,4) < le(i) & DisSav(:,4) >= le(i-1) ));
end
figure
AXM = [Dny(:,1),Dny(:,2)];
b = bar(le(2:end),AXM); 
b(1).FaceColor = 'k'; b(1).EdgeColor = 'k'; b(1).LineWidth = 2;
plot(ones(20,1)*mean(DisSav(:,1)),linspace(0,14,20),'-*r')
plot(linspace(mean(DisSav(:,1))-std(DisSav(:,1)),mean(DisSav(:,1))+std(DisSav(:,1)),20),ones(20,1)*40*mean(DisSav(:,1)),'--k')
b(2).FaceColor = 'r'; b(2).EdgeColor = 'r'; b(2).LineWidth = 2;
plot(ones(20,1)*mean(DisSav(:,2)),linspace(0,14,20),'-*r')
plot(linspace(mean(DisSav(:,2))-std(DisSav(:,2)),mean(DisSav(:,2))+std(DisSav(:,2)),20),ones(20,1)*40*mean(DisSav(:,2)),'--r')
goodplot;

figure
AXM = [Dny(:,3),Dny(:,2)];
b = bar(le(2:end),AXM); 
b(1).FaceColor = 'b'; b(1).EdgeColor = 'b'; b(1).LineWidth = 2;
plot(ones(20,1)*mean(DisSav(:,3)),linspace(0,14,20),'-*b')
plot(linspace(mean(DisSav(:,3))-std(DisSav(:,3)),mean(DisSav(:,3))+std(DisSav(:,3)),20),ones(20,1)*50*mean(DisSav(:,3)),'--b')
b(2).FaceColor = 'r'; b(2).EdgeColor = 'r'; b(2).LineWidth = 2;
plot(ones(20,1)*mean(DisSav(:,2)),linspace(0,14,20),'-*r')
plot(linspace(mean(DisSav(:,2))-std(DisSav(:,2)),mean(DisSav(:,2))+std(DisSav(:,2)),20),ones(20,1)*40*mean(DisSav(:,2)),'--r')
goodplot;



gg = hot(6);
figure; hold on
h = histogram(DisSav(:,1)); hold on
h.BinEdges = [0:0.06:2.5];
h.Normalization = 'PDF';
h.FaceColor = 'k';
h.EdgeColor = 'k';
set(h,'facealpha',0.1);
h.LineWidth = 2;
plot(ones(20,1)*mean(DisSav(:,1)),linspace(0,14,20),'-*k')
plot(linspace(mean(DisSav(:,1))-std(DisSav(:,1)),mean(DisSav(:,1))+std(DisSav(:,1)),20),ones(20,1)*6*mean(DisSav(:,1)),'--k')


h = histogram(DisSav(:,2)); hold on
h.BinEdges = [0:0.06:2.5];
h.Normalization = 'PDF';
h.FaceColor =  'r';
h.EdgeColor = 'r';
set(h,'facealpha',0.1);
h.LineWidth = 2;
plot(ones(20,1)*mean(DisSav(:,2)),linspace(0,14,20),'-*r')
plot(linspace(mean(DisSav(:,2))-std(DisSav(:,2)),mean(DisSav(:,2))+std(DisSav(:,2)),20),ones(20,1)*8*mean(DisSav(:,2)),'--r')
goodplot;

figure; hold on
h = histogram(DisSav(:,3)); hold on
h.BinEdges = [0:0.06:2.5];
h.Normalization = 'PDF';
h.FaceColor = 'b';
h.EdgeColor = 'b';
set(h,'facealpha',0.1);
h.LineWidth = 2;
plot(ones(20,1)*mean(DisSav(:,3)),linspace(0,14,20),'-*b')
plot(linspace(mean(DisSav(:,3))-std(DisSav(:,3)),mean(DisSav(:,3))+std(DisSav(:,3)),20),ones(20,1)*6*mean(DisSav(:,3)),'--b')

h = histogram(DisSav(:,2)); hold on
h.BinEdges = [0:0.06:2.5];
h.Normalization = 'PDF';
h.FaceColor =  'r';
h.EdgeColor =  'r';
set(h,'facealpha',0.1);
h.LineWidth = 2;
plot(ones(20,1)*mean(DisSav(:,2)),linspace(0,14,20),'-*r')
plot(linspace(mean(DisSav(:,2))-std(DisSav(:,2)),mean(DisSav(:,2))+std(DisSav(:,2)),20),ones(20,1)*8*mean(DisSav(:,2)),'--r')
goodplot;





