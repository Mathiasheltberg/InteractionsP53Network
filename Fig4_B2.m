clear all; close all; clc

%C1 = load('Sheng1.mat');
C1 = load('MdmxsiRNA_UV16.mat');

AlAf_1  = zeros(101,143);
cg1 = 0;
Peaks1 = zeros(30,2);
Before = zeros(50,101);
During = zeros(143,101);
After = zeros(143,101);
figure
for i = 1:15:100
    i
   %F1 = C1.UV8MX_data(1,i).YFP;
    F1 = C1.UV16MX_data(1,i).YFP;
    M = 50;
    if (F1(end) > 0)
        cg1 = cg1 + 1;
        g = find(F1(50:60) == max(F1(50:60)))+49;
        Peaks1(cg1,2) = F1(g); Peaks1(cg1,1) = g;
        Before(:,cg1) = F1(1:50)';
        L = find(F1 < mean(F1(1:50))); LL = L(min(find(L > g)));
        During(1:LL-M+1,cg1) = F1(M:LL);
        After(1:143-LL+1,cg1) = F1(LL:end);
        AlAf_1(cg1,1:length(F1)) = F1;
                
        subplot(1,2,cg1)
    l1 = linspace(1,50/3,50);
    plot(l1,Before(:,cg1),'LineWidth',3); hold on
    l2 = linspace(M/3,LL/3,LL-M+1);
    plot(l2,During(1:LL-M+1,cg1),'LineWidth',3);
    l3 = linspace(LL/3,142/3,143-LL);
    plot(l3,After(1:143-LL,cg1),'LineWidth',3);
    plot(ones(30,1)*50/3,linspace(0,250,30),'.k','MarkerSize',12)
    plot(linspace(0,143/3,30),ones(30,1)*mean(Before(:,cg1)),'--','color',[0.1 0.5 0.1])
    goodplot                
    end        
end
