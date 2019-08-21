clear all; close all; clc

P1 = importdata('Control_P.txt');
M1 = importdata('Control_M.txt');

P2 = importdata('XKD_P.txt');
M2 = importdata('XKD_M.txt');
H1 = M1./P1;
H2 = M2./P2;

figure;
le = 0.0:0.05:2.5;
Dny = zeros(length(le)-1,2);
Dny(1,1) = length(find(H1 < le(1) ) );
Dny(1,2) = length(find(H2 < le(1) ) );
for i = 2:length(le)-1;
    Dny(i,1) = length(find(H1 < le(i) & H1 >= le(i-1) ));
    Dny(i,2) = length(find(H2 < le(i) & H2 >= le(i-1) ));
end
AXM = [Dny(:,1), Dny(:,2)]; b = bar(le(2:end),AXM); b(1).FaceColor = 'k'; b(1).EdgeColor = 'k'; b(1).LineWidth = 2;
legend(['Control ',num2str(mean(H1)),' s ',num2str(std(H1))],['X-KD ',num2str(mean(H2)),' s ',num2str(std(H2))]); 
b(2).FaceColor = 'r'; b(2).EdgeColor = 'r'; b(2).LineWidth = 2;
goodplot;


figure
loglog(P1(1:500),M1(1:500),'.b','Markersize',15); hold on
loglog(P2(1:500),M2(1:500),'.r','Markersize',15); goodplot