clear all; close all; clc
C1 = load('Sheng1.mat');

AlAf_1  = zeros(101,143);
cg1 = 0;
Peaks1 = zeros(30,2);
Before = zeros(50,101);
During = zeros(143,101);
After = zeros(143,101);
CTot = 0; Amps = [];
for i = 1:101
    i
   F1 = C1.UV8MX_data(1,i).YFP;
    M = 50;
    
    if (F1(end) > 0)
        cg1 = cg1 + 1

        As = []; Ps = []; Bs = []; Cs = 0; Ct = 0;
        for u = 51:length(F1)-1
            if (F1(u) > F1(u-1) && F1(u) > F1(u+1))
                if (u <= 3)
                    if (F1(u) == max(F1(1:u+3)))
                        Cs = Cs+1; As(Cs) = F1(u); Ps(Cs) = u;
                    end
                elseif (u >= 141)
                    if (F1(u) == max(F1(u-3:end)))
                        Cs = Cs+1; As(Cs) = F1(u); Ps(Cs) = u;
                    end
                else
                    if (F1(u) == max(F1(u-3:u+3)))
                        Cs = Cs+1; As(Cs) = F1(u); Ps(Cs) = u;
                    end
                end
            end
        end
        Ch = find(Ps(2:end)-Ps(1:end-1) <= 5);
        As(Ch) - As(Ch+1);
        for u = 1:length(As)-1;
            CTot = CTot +1;
            if u == 1               
                R = Ps(u):Ps(u+1);
                AU = find(F1(R) == min(F1(R)));                
                Amps(CTot) = As(u)-F1(AU+Ps(u));
            else
                R1 = Ps(u):Ps(u+1);
                AU1 = find(F1(R1) == min(F1(R1))); 
                R2 = Ps(u-1):Ps(u);
                AU2 = find(F1(R2) == min(F1(R2))); 
                Amps(CTot) = As(u)- mean([F1(AU1+Ps(u)) F1(AU2+Ps(u-1))]);
            end
        end
        
        g = find(F1(50:end) == max(F1(50:end)))+49;
        Peaks1(cg1,2) = F1(g); Peaks1(cg1,1) = g;
        Before(:,cg1) = F1(1:50)';
        L = find(F1 < mean(F1)); LL = L(min(find(L > g)));
        During(1:LL-M+1,cg1) = F1(M:LL);
        Nul = min(find(F1 == 0));
        if (length(Nul) == 0)
            Nul = 143;
        end
        After(1:Nul-50,cg1) = F1(51:Nul);
        
        AlAf_1(cg1,1:length(F1)) = F1;
    end
end

figure
B = sort(Amps);
C = min(find(B > 120));
h = histfit(B(2:C)/2,17,'Gamma'); goodplot
set(h(1),'facecolor',[0.9 0.9 0.99]); 
set(h(1),'edgecolor',[0.1 0.1 0.5],'LineWidth',4); 
set(h(2),'color',[0.9 0.2 0.0],'LineWidth',4)