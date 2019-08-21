clear all; close all; clc

C1 = load('Sheng1.mat');
C2 = load('Sheng2.mat');



AlAf_1  = zeros(101,143);
AlAf_2  = zeros(101,143);
cg1 = 0; cg2 = 0;
Peaks1 = zeros(101,2);
Peaks2 = zeros(101,2);
Before = zeros(50,101);
During = zeros(143,101);
After = zeros(143,101);
l = linspace(1/3,143/3,143);
for i = 1:101
    i
    F1 = C1.UV8MX_data(1,i).YFP;
    F2 = C2.UV8cnt_data(1,i).YFP;
    
    M = 50;
    if (F1(end) > 0)
                cg1 = cg1 + 1;
         plot(l,F1,'.b'); hold on; 
         AlAf_1(cg1,1:length(F1)) = F1; 

        g = find(F1(50:end) == max(F1(50:end)))+49;
        Peaks1(cg1,2) = F1(g); Peaks1(cg1,1) = (F1(g)-F1(50))/(g-50);
    end
    if (F2(end) > 0)
        cg2 = cg2 + 1;
        AlAf_2(cg2,1:length(F2)) = F2; 
        plot(l,F2,'.r')
        g = find(F2(50:end) == max(F2(50:end)))+49;
        Peaks2(cg2,2) = F2(g); Peaks2(cg2,1) = (F2(59)-F2(50))/(59-50);
    end
        
end

Means1 = zeros(143,1);
Means2 = zeros(143,1);
for i = 1:143
    ML1 = find(AlAf_1(:,i));
    ML2 = find(AlAf_2(:,i));
    Means1(i,1) = mean(AlAf_1(ML1,i));
    Means2(i,1) = mean(AlAf_2(ML2,i));
end


plot(l,Means1,'LineWidth',5); goodplot
plot(l,Means2,'LineWidth',5); goodplot
axis([0 50 0 500])


