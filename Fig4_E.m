clear all; close all; clc

C1 = load('Sheng1.mat');
C2 = load('Sheng2.mat');

AlAf_1  = zeros(101,143);
cg1 = 0; cg2 = 0;
Peaks1 = zeros(101,2);
Peaks2 = zeros(101,2);
Before = zeros(50,101);
During = zeros(143,101);
After = zeros(143,101);

for i = 1:3
    i
    F1 = C1.UV8MX_data(1,i).YFP;
    F2 = C2.UV8cnt_data(1,i).YFP;
    M = 50; DD = 5;
    if (F1(end) > 0)
        plot(F1,'b'); hold on; 
        cg1 = cg1 + 1;
        g = find(F1(50:80) == max(F1(50:80)))+49;
        Peaks1(cg1,2) = F1(g); Peaks1(cg1,1) = (F1(g)-F1(g-D))/(D);
        plot(g,F1(g),'*b')
        plot(g-D,F1(g-D),'ob')
    end
    if (F2(end) > 0)
        plot(F2,'r'); hold on
        cg2 = cg2 + 1;
        g = find(F2(50:80) == max(F2(50:80)))+49;
        Peaks2(cg2,2) = F2(g); Peaks2(cg2,1) = (F2(g)-F2(g-D))/(D);
        plot(g,F2(g),'*r')
        plot(g-3,F2(g-3),'or')
    end
        
end

figure
B1 = mean(Peaks1(1:cg1,1)); eB1 = std(Peaks1(1:cg1,1))/sqrt(cg1);
B2 = mean(Peaks2(1:cg2,1)); eB2 = std(Peaks2(1:cg2,1))/sqrt(cg2);
b = bar(1:2,[B1 B2]); hold on
b(1).FaceColor = [0.0 0.7 0.];
b(1).EdgeColor = [0.1 0.5 0.1];
b.LineWidth = 4;
be = errorbar(1:2,[B1 B2],[(eB1) (eB2)],'*');
be.LineWidth = 4;
be.Color = [0.1 0.5 0.1];
goodplot
