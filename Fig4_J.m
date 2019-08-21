clear all; close all; clc

C1 = load('MdmxsiRNA_UV16.mat');

AlAf_1  = zeros(101,143);
cg1 = 0;
Peaks1 = zeros(30,2);
Before = zeros(50,101);
During = zeros(143,101);
After = zeros(143,101);

for i = 1:101
    i
    F1 = C1.UV16MX_data(1,i).YFP;
    M = 50;
    if (F1(end) > 0)
        cg1 = cg1 + 1;
        g = find(F1(50:end) == max(F1(50:end)))+49;
        Peaks1(cg1,2) = F1(g); Peaks1(cg1,1) = g;
        Before(:,cg1) = F1(1:50)';
        L = find(F1 < mean(F1)); LL = L(min(find(L > g)));
        During(1:LL-M+1,cg1) = F1(M:LL);
        After(1:143-LL+1,cg1) = F1(LL:end);
        AlAf_1(cg1,1:length(F1)) = F1;
        
        if (LL < 95)
            Fs = 3;  T = 1/Fs;
            X1 = After(1:143-LL+1,cg1);
            Y1 = fft(X1)/1;
            L = length(X1);
            P2 = (abs(Y1/L)); P1 = P2(1:L/2+1);
            P1(2:end-1) = 2*P1(2:end-1);
            f = Fs*(0:(L/2))/L; [apD bpD] = max(P1(3:end));
            plot(f,P1,'.','color',[0.1 0.1 0.5],'MarkerSize',15); hold on; plot(f(bpD+2),apD,'ro');
            axis([0.05 1 0 40])
            goodplot
        end
        %
        %     X1 = Before(1:50,cg1);
        %     Y1 = fft(X1)/1;
        %     L = length(X1);
        %     P2 = (abs(Y1/L)); P1 = P2(1:L/2+1);
        %     P1(2:end-1) = 2*P1(2:end-1);
        %     f = Fs*(0:(L/2))/L; [apD bpD] = max(P1(3:end));
        %     plot(f,P1,'.','color',[0.1 0.1 0.5],'MarkerSize',15); hold on; plot(f(bpD+2),apD,'ro');
        %     goodplot
        %
%         figure(2)
%         subplot(1,2,cg1)
%     l1 = linspace(1,50/3,50);
%     plot(l1,Before(:,cg1),'LineWidth',3); hold on
%     l2 = linspace(M/3,LL/3,LL-M+1);
%     plot(l2,During(1:LL-M+1,cg1),'LineWidth',3);
%     l3 = linspace(LL/3,142/3,143-LL);
%     plot(l3,After(1:143-LL,cg1),'LineWidth',3);
%     plot(ones(30,1)*50/3,linspace(0,250,30),'.k','MarkerSize',12)
%     plot(linspace(0,143/3,30),ones(30,1)*mean(Before(:,cg1)),'--','color',[0.1 0.5 0.1])
%     goodplot
        
        
    end
    

    
end

% figure;
% 
% % plot(AlAf_1(:,:)','-.'); hold on;
% l1 = linspace(1,50/3,50);
% plot(l1,Before(:,1),'LineWidth',3); hold on
% l2 = linspace(M/3,LL/3,LL-M+1);
% plot(l2,During(1:LL-M+1,1),'LineWidth',3);
% l3 = linspace(LL/3,142/3,143-LL);
% plot(l3,After(1:143-LL,1),'LineWidth',3);
% plot(ones(30,1)*50/3,linspace(0,250,30),'.k','MarkerSize',12)
% plot(linspace(0,143/3,30),ones(30,1)*mean(Before(:,1)),'--','color',[0.1 0.5 0.1])
% goodplot