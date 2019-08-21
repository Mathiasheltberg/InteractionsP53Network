clear all; close all; clc
%%% Here we calculate the Fourier Spectrum of oscillations
A = load('MdmxsiRNA_data.mat');
Before = zeros(170,96);
During = zeros(170,96);
After = zeros(170,96);
cg = 0;
Cel = zeros(170*2+1,96);
CH = zeros(96,1);
x = linspace(-170,170,341);
AllP = [];
for i = 1:101
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
        
        si = 171-M+1; sf = 170*2-M+1;
        Before(end-M+1:end,cg) = F1(1:M);
        L = find(F1 < mean(F1)); LL = L(min(find(L > g)));
        During(1:LL-M,cg) = F1(M:LL-1);
        After(1:170-LL,cg) = F1(LL+1:end);
        Cel(si:sf,cg) = F1';
        CH(cg,1) = max(During(:,cg));
        
        l = linspace(0,length(After(1:170-LL))/3,length(After(1:170-LL)));
        l2 = linspace(0,length(After(1:170-LL))/3,10*length(After(1:170-LL)));
        
        
        Fs = 3;  T = 1/Fs;
        X1 = After(1:170-LL,1);
        Y1 = fft(X1)/1;
        L = length(X1);
        P2 = (abs(Y1/L)); P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs*(0:(L/2))/L; [apD bpD] = max(P1(10:end));
        plot(f,P1,'.','color',[0.1 0.5 0.1],'MarkerSize',15); hold on; plot(f(bpD+9),apD,'ro');
        AllP = [AllP;1/f(bpD+9)];
        goodplot
    end
    
end
goodplot;
figure
h = histfit(AllP,10,'Gamma'); goodplot
set(h(1),'facecolor',[0.9 0.99 0.9]); 
set(h(1),'edgecolor',[0.1 0.5 0.1],'--','LineWidth',4); 
set(h(2),'color',[0.1 0.9 0.0],'LineWidth',4)

