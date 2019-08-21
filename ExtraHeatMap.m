clear all; close all; clc
%%% Here we calculate the Fourier Spectrum of oscillations
A = load('MdmxsiRNA_data.mat');
Before = zeros(170,96);
During = zeros(170,96);
After = zeros(170,96);
All = zeros(120,96);
cg = 0;
Cel = zeros(170*2+1,96);
CH = zeros(96,1);
x = linspace(-170,170,341);
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
        xg = find(F1 > 150);
        if (xg > 0);
            F1(xg) = 150;
        end
%         F1 = F1/mean(F1);
        All(:,cg) = F1(1:120);
        plot(F1); hold on
    end
    
end
plot(mean(All'),'LineWidth',3)

figure
imagesc(All)
colorbar
colormap('jet')
goodplot