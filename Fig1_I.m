clear all; close all; clc
%%% Here we calculate the Fourier Spectrum of oscillations
A = load('MdmxsiRNA_data.mat');
Before = zeros(170,96);
During = zeros(170,96);
After = zeros(170,96);
All = zeros(170,96);
cg = 0;
Cel = zeros(170*2+1,96);
CH = zeros(96,1);
x = linspace(-170,170,341);
for i = 76:76
%for i = 1:101
    i;
    F1m = A.data(1,i).mitosis{1};
    F1 = A.data(:,i).YFP;
    if (F1(end)>0)
        M = min(F1m);
        if (M == 0)
            M = 1;
        end
        
        cg = cg+1;
        g = find(F1(F1m(1)+1:F1m(1)+15) == max(F1(F1m(1)+1:F1m(1)+15)))+F1m(1)+1;
        All(:,cg) = F1;
        si = 171-M+1; sf = 170*2-M+1;
        Before(end-M+1:end,cg) = F1(1:M);
        L = find(F1 < mean(F1)); LL = L(min(find(L > g)));
        During(1:LL-M,cg) = F1(M:LL-1);
        After(1:170-LL,cg) = F1(LL+1:end);
        Cel(si:sf,cg) = F1';
        CH(cg,1) = max(During(:,cg));
        figure
        l = linspace(0,length(After(1:170-LL))/3,length(After(1:170-LL)));
        l2 = linspace(0,length(After(1:170-LL))/3,10*length(After(1:170-LL)));
        plot(l,After(1:170-LL,cg),'color',[0.1 0.5 0.1],'LineWidth',5); hold on
       % plot(l2,20 + 0.03*(l2-20).^2 + 10*sin(2*pi*0.3125*l2+5.0),'-','LineWidth',3)
        plot(l2,25 + 10*sin(2*pi*0.3125*l2+5.0),'-','LineWidth',3)
        goodplot

    end
    
end
goodplot;

figure
imagesc(All)
