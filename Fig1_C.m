clear all; close all; clc

A = load('MdmxsiRNA_data.mat');
Before = zeros(170,96);
During = zeros(170,96);
After = zeros(170,96);
cg = 0;
Cel = zeros(170*2+1,96);
CM = zeros(170*2+1,1);
x = linspace(-170/3,170/3,341);

for i = 1:101
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
        si = 171-M; sf = 170*2-M;
        L = find(F1 < mean(F1)); LL = L(min(find(L > g)))-1;
        Cel(si:sf,cg) = F1';
        plot(x,Cel(:,cg),'.'); hold on;
    end
    
end

for j = 1:170*2+1
    L = find(Cel(j,:) > 0);
    if (length(L) > 0)
        CM(j) = mean(Cel(j,L));
    end
end
plot(x,CM,'color',[0.1 0.5 0.1],'LineWidth',5); goodplot
goodplot
