clear all; close all; clc
Val = zeros(100,100);
for test1 = 1:100;
    for test2 = 1:100;
        l2 = test2/20.; l1 = test1/10.;    
        k1 = 0.1; 
  
        k2 = 1.0; k4 = 0.01;
        k7 = 0.13; k8 = 0.1;

        k2 = k2*(1 + l1*1);
        k4 = k4*(1 + l2*1);
        K = k1*k8/(k7*k2);
        fpp = K/2.0*(1+sqrt(1+4*k4/K));
        Val(test1,test2) =  0.0864/fpp;
        fpm = k7/k8*fpp;

       
    
    end
end

figure
imagesc(Val);
colormap('hot'); colorbar
goodplot