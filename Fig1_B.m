clear all; close all; clc
%%% Load data
A = load('MdmxsiRNA_data.mat');


Before = zeros(101,1); %Before will be the data before Mdmx suppression
After = zeros(101,1); %After will be the data after Mdmx suppression
AlAf  = zeros(101,170); %This will be the data for all cells
cg = 0;
R = [40 59 64 95];
for i = 1:101
    i;
    F1m = A.data(1,i).mitosis{1};
    F1 = A.data(:,i).YFP;   
    
    if (F1(end)>0 && length(find(R == i))>0)
        x = linspace(-F1m(1)/2,(170-F1m(1)-1)/2,170);     %% This aligns the cells    
        cg = cg+1;  
        subplot(2,2,cg)
        g = find(F1(F1m(1)+1:F1m(1)+15) == max(F1(F1m(1)+1:F1m(1)+15)))+F1m(1)+1;        
        plot(x,F1,'color',[0.1 0.5 0.1],'LineWidth',3); hold on; plot(x(g-1),F1(g-1),'k*'); hold on
        plot(zeros(40,1),linspace(0,200,40),'.k'); axis([-20 120 0 200]);
        axis([-10 80 0 175])
        goodplot
    end
    
end


