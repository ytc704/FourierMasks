% Yutong Chen, Princeton University
% Afonso Bandeira, Princeton University
% Dustin Mixon, Air Force Institute of Technology
%
% Free to use. Please cite our paper "Phase retrieval from power spectra of
% masked signals" if the code is used in publications.

% demo
clear;
clc;

%% signal
M = 64; % signal length
x = random('Normal',0,1,M,1);

%% polarization algorithm

k = 3; % number of masks
VarNoiseM2 = 0.1; % Gaussian noise variance
fracRmv = 0.01; % fraction of edges to be removed from graph

[t_Pol,relErr_Pol,newMu,xhat_Pol] = testPolarization(x,k,VarNoiseM2,fracRmv);


%% plot

subplot(2,1,1)
plot(1:M,x);
axis([0 M+1 min(x)-1 max(x)+1])
title('True signal');

subplot(2,1,2)
plot(1:M,real(xhat_Pol));
axis([0 M+1 min(x)-1 max(x)+1])
title({'Reconstruction via polarization',['Time: ',num2str(t_Pol),'s    L2 Error: ',num2str(100*relErr_Pol),'%']})