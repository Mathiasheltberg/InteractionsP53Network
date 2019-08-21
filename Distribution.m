clear all; close all; clc
MM = 300;
DisSav = zeros(MM,2);
Tmax = 300; ts = 0.1; dt = 0.001; N = 10;
Pe = zeros(MM,1); PIAL = zeros(round(Tmax/ts),MM);
enois = 0.06; inois = 0.06;
for testN = 1:MM
    testN
    NN = round(Tmax/dt);
    Sav = zeros(round(Tmax/ts),3);
    P_I = ones(NN,1)*0.0322;       M = ones(NN,1)*0.0349;
    d1 = 0; d2 = 0; d3 = 0; d4 = 0; d5 = 0; d6 = 0;
    k1 = 0.1; TDel1 = round(2000);
    k2 = 1.0; k4 = 0.01;
    k7 = 0.13; k8 = 0.1;
    k10 = k1; k20 = k2; k80 = k8; k70 = k7; k40 = k4;
    TDel10 = TDel1;
    
    RT = 0; click = 1; RT0 = 0; c = 1;
    
    X = 1.;
    lambda1 = 3.; lambda2 = 0.; lambda3 = 0.15;
    
    k20 = k20*(1 + lambda1*X);
    k40 = k40*(1 + lambda2*X);
    k70 = k70*(1 - lambda3*X);
    if (testN < MM)
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
        
        if (testN < MM)
            P_I(c) = P_I(c-1) + dt*(dPi) + sqrt(dt)*(sqrt(k1) - sqrt(k2*M(c-1)*P_I(c-1)/(k4 + P_I(c-1))))*inois*randn;
            M(c) = M(c-1) + dt*(dM) + sqrt(dt)*(sqrt( k7*(P_I(mod(c-TDel1,NN)+1))) - sqrt(k8*M(c-1)))*inois*randn;
        else
            P_I(c) = P_I(c-1) + dt*(dPi); 
            M(c) = M(c-1) + dt*(dM);
        end
        
        %%%%%%%%%%%%% Make sudden impace
        if (RT > 50 && RT0 < 50);
            DisSav(testN,1) = M(c)/P_I(c);
            k2 = k20/(1 + lambda1*X); k4 = k40/(1 + lambda2*X); k7 = k70/(1 - lambda3*X);
            click0 = click; end
        
        %%%%%%%%%%%% Secure Against Runaways
        if (P_I(c) < 0); P_I(c) = 0.000001; end;
        if (M(c) < 0);   M(c) = 0.000001; end;
        
        %%%%%%%%%%%% Save elements
        if (RT > ts*click);
            Sav(click,1) = RT;
            Sav(click,2) = P_I(c);
            Sav(click,3) = M(c);
            PIAL(click,testN) = P_I(c);
            click=click+1; end
        
        if (RT > 150 && RT0 < 150);
            DisSav(testN,3) = M(c)/P_I(c); end
        
        if (RT > 100 && RT0 < 100);
            DisSav(testN,2) = M(c)/P_I(c); end
        RT0 = RT;
        
    end
    if (testN < MM)
        plot(Sav(1:end-1,1),Sav(1:end-1,2),'.'); hold on;
    end
    DisSav(testN,4) = M(c)/P_I(c);
    Pe(testN,1) = max(P_I);
end
YC = mean(PIAL');
plot(Sav(1:end-1,1),Sav(1:end-1,2),'r','LineWidth',3); hold on;
%plot(Sav(1:end-1,1),YC(1:end-1),'LineWidth',3); hold on;
goodplot
axis([0.0 Sav(end-1,1) 0.0 max(Sav(:,2))*2])
goodplot



XX = [DisSav(:,2); DisSav(:,3); DisSav(:,4)];

figure;
le = 0.0:0.1:2.5;
Dny = zeros(length(le)-1,2);
Dny(1,1) = length(find(DisSav(:,1) < le(1) ) );
Dny(1,2) = length(find(DisSav(:,4) < le(1) ) );
for i = 2:length(le)-1;
    Dny(i,1) = length(find(DisSav(:,1) < le(i) & DisSav(:,1) >= le(i-1) ));
    Dny(i,2) = length(find(DisSav(:,4) < le(i) & DisSav(:,4) >= le(i-1) ));
end
AXM = [Dny(:,1), Dny(:,2)]; b = bar(le(2:end),AXM); b(1).FaceColor = 'k'; b(1).EdgeColor = 'k'; b(1).LineWidth = 2;
%legend(['Control ',num2str(mean(DisSav(:,1))),' s ',num2str(std(DisSav(:,1)))],['X-KD ',num2str(mean(DisSav(:,2))),' s ',num2str(std(DisSav(:,2)))]);
b(2).FaceColor = 'r'; b(2).EdgeColor = 'r'; b(2).LineWidth = 2;
goodplot;

figure
h = histfit(Pe,10,'Gamma'); goodplot
% set(h(1),'facecolor',[0.9 0.99 0.9]);
% set(h(1),'edgecolor',[0.1 0.5 0.1],'--','LineWidth',4);
% set(h(2),'color',[0.1 0.9 0.0],'LineWidth',4)

figure
h = histogram(XX);
h.BinEdges = [0:0.1:3];
h.Normalization = 'PDF';
h.FaceColor = [0.0 0.0 0.7];
h.EdgeColor = [0.1 0.1 0.5];
h.LineWidth = 4;
hold on
h = histogram(DisSav(:,1));
h.BinEdges = [0:0.1:3];
h.Normalization = 'PDF';
h.FaceColor = [0.7 0.0 0.0];
h.EdgeColor = [0.5 0.1 0.5];
h.LineWidth = 4;
goodplot



