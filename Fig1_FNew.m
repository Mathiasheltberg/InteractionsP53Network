clear all; close all; clc
%%%%%% To Calculate Height of Pulse
A = load('MdmxsiRNA_data.mat');
Before = zeros(170,96);
During = zeros(170,96);
After = zeros(170,96);
cg = 0;
Cel = zeros(170*2+1,96);
CH = zeros(96,1);
SL = zeros(96,1);
SL2 = zeros(96,1);

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
        
        si = 171-M+1; sf = 170*2-M+1;
        Before(end-M+1:end,cg) = F1(1:M);
        L = find(F1 < mean(F1)); LL = L(min(find(L > g)));
        During(1:LL-M,cg) = F1(M:LL-1);
        After(1:170-LL,cg) = F1(LL+1:end);
        Cel(si:sf,cg) = F1';
%         figure
        a1 = find(F1 == max(During(:,cg)))
        a2 = LL;
%         
        CH(cg,1) = max(During(:,cg));
        SL(cg,1) = (F1(a1)-F1(a1-3))/3.;
        SL2(cg,1) = (F1(a1)-F1(a1+3))/3.;
        if (SL(cg,1) < 0)
            figure; plot(F1); hold on; plot(M,F1(M),'*');  plot(a2,F1(a2),'*'); 
        end
        end    
end

figure; hold on
B = sort(CH);
[p,ci] = gamfit(CH);
h = histogram(CH);
h.BinEdges = [0:25:250];
h.Normalization = 'count';
h.EdgeColor = [0.0 0.7 0.0];
h.FaceColor = [0.6 0.9 0.7];
h.LineWidth = 4;
goodplot
y = gampdf(linspace(0,250,1001),p(1),p(2));
plot(linspace(0,250,1001),y,'-','color',[0.8 0.0 0.3],'LineWidth',5)

figure; hold on
B = sort(SL);
[p,ci] = gamfit(SL);
h = histogram(SL);
h.BinEdges = [0:3.0:40];
h.Normalization = 'number';
h.EdgeColor = [0.0 0.5 0.0];
h.FaceColor = [0.2 0.6 0.2];
h.LineWidth = 4;
goodplot
y = gampdf(linspace(0,45,1001),p(1),p(2));
plot(linspace(0,45,1001),y,'-','color',[0.8 0.0 0.3],'LineWidth',5)

B1 = mean(SL); eB1 = std(SL)/sqrt(length(SL));
B2 = mean(SL2); eB2 = std(SL2)/sqrt(length(SL2));

figure
b = bar(1:2,[B1 B2]); hold on
b(1).FaceColor = [0.0 0.7 0.0];
b(1).EdgeColor = [0.6 0.9 0.7];
b.LineWidth = 4;
be = errorbar(1:2,[B1 B2],[eB1 eB2],'*');
be.LineWidth = 4;
be.Color = [0.1 0.5 0.1];
goodplot
