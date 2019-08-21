clear all; close all; clc
%%%% This is to find the mean slope
A = load('MdmxsiRNA_data.mat');
Before = zeros(170,96);
During = zeros(170,96);
After = zeros(170,96);
cg = 0;
Cel = zeros(170*2+1,96);
CM = zeros(170*2+1,1);
xl = linspace(-170,170,341);
cg1 = 0; cg2 = 0;
Peaks1 = zeros(101,2);
for i = 1:101
    i;
    F1m = A.data(1,i).mitosis{1};
    F1 = A.data(:,i).YFP;
    if (F1(end)>0)
        M = min(F1m)+1;
        if (M == 0)
            M = 1;
        end
       %figure
        cg = cg+1;
        g = find(F1(F1m(1)+1:F1m(1)+15) == max(F1(F1m(1)+1:F1m(1)+15)))+F1m(1)+1;
        plot(F1); hold on;
        GG = M+1;
        plot(g-1,F1(g-1),'*');
        plot(M,F1(M),'o');
        si = 171-M; sf = 170*2-M;
        
        L = find(F1 < mean(F1)); LL = L(min(find(L > g)))-1;
        
        Cel(si:sf,cg) = F1';
        cg1 = cg1+1;
        
        Peaks1(cg1,1) = (F1(g)-F1(GG))/(g-1-(GG));
        Peaks2(cg1,1) = (F1(g)-F1(LL+1))/(g-1-LL+1);
        % plot(x,Cel(:,cg),'.'); hold on;
        
    end
    
end

for j = 1:170*2+1
    L = find(Cel(j,:) > 0);
    if (length(L) > 0)
        CM(j) = mean(Cel(j,L));
    end
    
    
end
% xl = linspace(-170/2,170/2,341);
% plot(xl,CM,'color',[0.1 0.5 0.1],'LineWidth',7); hold on; goodplot
% plot(zeros(100,1),linspace(0,120,100),'--k')
% axis([-10/3 20/3 20 100])
% [x,y] = ginput(3)
% plot(x(1:2),y(1:2),'color',[0.8 0.5 0.1],'LineWidth',4)
% plot(x(2:3),y(2:3),'color',[0.9 0.6 0.2],'LineWidth',4)
% figure
% TT = [(y(2)-y(1))/(x(2)-x(1)); (abs(y(3)-y(2)))/(x(3)-x(2))];
% bar(TT,'FaceColor',[0.9 .99 .9],'EdgeColor',[0.1 .5 .1],'LineWidth',5.5)
% axis([0 3 0 30])
% goodplot

figure
B1 = mean(Peaks1(1:cg1,1)); eB1 = std(Peaks1(1:cg1,1))/sqrt(cg1);
b = bar(1,[B1]); hold on
b(1).FaceColor = [0.0 0.7 0.];
b(1).EdgeColor = [0.1 0.5 0.1];
b.LineWidth = 4;
be = errorbar(1,[B1],[(eB1)],'*');
be.LineWidth = 4;
be.Color = [0.1 0.5 0.1];
goodplot
