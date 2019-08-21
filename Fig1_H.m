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
CTot = 0; Amps = [];
AllP = [];
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
        CH(cg,1) = max(During(:,cg));
        FN = F1(LL+1:end);
        
        plot(F1,'-dr'); hold on
        As = []; Ps = []; Bs = []; Cs = 0; Ct = 0;
        Per = [];
        for u = 2:length(FN)-1
            if (FN(u) > FN(u-1) && FN(u) > FN(u+1))
                if (u <= 3)
                    if (FN(u) == max(FN(1:u+3)))
                        Cs = Cs+1; As(Cs) = FN(u); Ps(Cs) = u;
                        Per = [Per;u];
                    end
                elseif (u >= length(FN)-3)
                    if (FN(u) == max(FN(u-3:end)))
                        Cs = Cs+1; As(Cs) = FN(u); Ps(Cs) = u;
                        Per = [Per;u];
                    end
                else
                    if (FN(u) == max(FN(u-3:u+3)))
                        Cs = Cs+1; As(Cs) = FN(u); Ps(Cs) = u;
                        Per = [Per;u];
                    end
                end
            end
        end
        AllP = [AllP; mean(Per(2:end)-Per(1:end-1))];
        Ch = find(Ps(2:end)-Ps(1:end-1) <= 5);
        As(Ch) - As(Ch+1);
        for u = 1:length(As)-1;
            CTot = CTot +1;
            if u == 1
                R = Ps(u):Ps(u+1);
                AU = find(FN(R) == min(FN(R)));
                Amps(CTot) = As(u)-FN(AU(1)+Ps(u));
            else
                R1 = Ps(u):Ps(u+1);
                AU1 = find(FN(R1) == min(FN(R1)));
                R2 = Ps(u-1):Ps(u);
                AU2 = find(FN(R2) == min(FN(R2)));
                Amps(CTot) = As(u)- mean([FN(AU1(1)+Ps(u)) FN(AU2(1)+Ps(u-1))]);
            end
        end
        
        
        
    end
    
end
goodplot;
AAA = sort(Amps/2);
figure
h = histfit(AAA(2:end-2),20,'Gamma'); goodplot
set(h(1),'facecolor',[0.9 0.99 0.9]); 
set(h(1),'edgecolor',[0.1 0.5 0.1],'LineWidth',4); 
set(h(2),'color',[0.9 0.2 0.0],'LineWidth',4)
pd = fitdist(AAA(2:end)','Gamma');
[h p] = chi2gof(AAA(2:end)','CDF',pd)
figure
h = histogram(Amps/2);
h.NumBins = 15;
h.BinEdges = [0:2:40];
h.Normalization = 'count';
h.FaceColor = [0.0 0.7 0.0];
h.EdgeColor = [0.1 0.5 0.1];
h.LineWidth = 4;
goodplot
[phat,pci] = gamfit(AAA(2:end-2))
mean(AAA(2:end-2))
