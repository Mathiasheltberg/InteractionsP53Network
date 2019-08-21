clear all; close all; clc
A = load('MdmxsiRNA_data.mat');
Af = []; cg = 0;
N = 101;
for i = 1:N
    i;
    F1m = A.data(1,i).mitosis{1};
    F1 = A.data(:,i).YFP;
    if (F1(end)>0)
        M = min(F1m);
        if (M == 0)
            M = 1;
        end
        cg = cg+1
        g = find(F1(F1m(1)+1:F1m(1)+15) == max(F1(F1m(1)+1:F1m(1)+15)))+F1m(1)+1;
        si = 171-M+1; sf = 170*2-M+1;
        L = find(F1 < mean(F1)); LL = L(min(find(L > g)));
        Af = [Af; {F1(LL+1:end)}];
    end
end
%%%%%% End import

figure
xi = 0;
for i = 1:length(Af)
    xo = xi + length(Af{i}); xl = linspace(xi,xo,xo-xi);
    plot(xl,Af{i}/mean(Af{i})); hold on; goodplot
    xi = xo;
end
%%%% Ends initial figure
%
figure
xxi = 1; DMed = []; xi = 0; NMed = [];
Firsts = [];
gg = jet(N)
for i = 1:length(Af)
    i
    xo = xi + length(Af{i}); xl = linspace(xi,xo,xo-xi);
    plot(xl,Af{i}./mean(Af{i}),'LineWidth',1); hold on
    xi = xo;
    X = Af{i}/mean(Af{i});
    As = []; Ps = []; Bs = []; Cs = 0; Ct = 0;
    for u = 2:length(X)-1
        if (X(u) > X(u-1) && X(u) > X(u+1))
            if (u <= 3)
                if (X(u) == max(X(1:u+3)))
                    Cs = Cs+1; As(Cs) = X(u); Ps(Cs) = u;
                end
            elseif (u >= length(X)-3)
                if (X(u) == max(X(u-3:end)))
                    Cs = Cs+1; As(Cs) = X(u); Ps(Cs) = u;
                end
            else
                if (X(u) == max(X(u-3:u+3)))
                    Cs = Cs+1; As(Cs) = X(u); Ps(Cs) = u;
                end
            end
        end
    end
    
    plot(xl(Ps),As,'*')
    X = X(Ps(1):Ps(end));
    
    xxo = xxi + length(X);
    xxl = linspace(xxi,xxo,xxo-xxi);
    xxi = xxo;
    

   % plot(xxl,X,'-'); hold on
    
    if (i>1)
        DMed(end) = mean([DMed(end),X(1)]);
        X = X(2:end);        
    end
    DMed = [DMed; X];
end
goodplot
figure
plot(DMed); hold on
goodplot;

Med = []; nd = 6; ndf = 6; All = [];
for u = 1:DMed
    for i = nd+1:length(DMed)-ndf
        Med = [Med; mean(DMed(i-nd:i+ndf))];
    end
    All = [All; DMed(nd+1:end-ndf)];
end
figure
plot(All); hold on
plot(Med); goodplot
%%%% Ends calculating runnign mean

% NMed = []; nd = 2; ndf = 2; NAll = [];
% 
% for i = nd+1:length(All)-ndf
%     NMed = [NMed; mean(All(i-nd:i+ndf))];
% end
% NAll = [NAll; All(nd+1:end-ndf)];
% figure
% plot(NAll,'g'); hold on; plot(NMed,'r'); goodplot

figure
TP = 0; CC =[];
y = (All-Med);
% y = y/max(y);
plot(y); hold on; plot(Med);goodplot
for i = 1:15
    TP = TP + 1;
    y1 = y(1:end-TP+1);
    y2 = y(TP:end);
    CC = [CC;abs(corr(y1,y2))];
end

D1 = min(find(CC < 0.5)); TP = D1;
TP = TP + 1;
y1 = y(1:end-TP+1);
y2 = y(TP:end);
%

figure
plot(log(y1),log(y2),'.b'); hold on
plot(log(y1),log(y2),'b'); hold on
goodplot

figure
plot((y1),(y2),'.b'); hold on
plot((y1),(y2),'b'); hold on
goodplot
% %

%
z1 = (y1);
z2 = (y2);
N = 50;
zz1 = linspace(-1,1,N+1);
zz2 = linspace(-1,1,N+1);
AA = zeros(N,N);
for i = 2:N+1
    for j = 2:N+1
        n = find(z1 > zz1(i-1) & z1 < zz1(i) & z2 > zz2(j-1) & z2 < zz1(j));
        if (length(n) == 0)
            AA(i-1,j-1) = 0;
        else
            AA(i-1,j-1) = length(n);
%             if (length(n) > 10) 
%                 AA(i-1,j-1) = 10;
%             end
        end
    end
end
figure
imagesc(AA); goodplot
colorbar
colormap('jet')

