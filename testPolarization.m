function [t,RelativeError,newMu,xhat] = testPolarization(x,k,VarNoiseM2,fracRmv)
% testPolarization See 'makeMeasurements' and 'reconstructPolarized'.
%
% Yutong Chen, Princeton University
% Afonso Bandeira, Princeton University
% Dustin Mixon, Air Force Institute of Technology
%
% Free to use. Please cite our paper "Phase retrieval from power spectra of
% masked signals" if the code is used in publications.

VarNoise = VarNoiseM2/(length(x)^2); % noise level

%% make masks and measurements

[W,Wphase,y,phi,spg,newMu] = makeMeasurements(x,k,VarNoise);

fprintf('\n');
disp('========POLARIZATION========');
fprintf('\n');
disp(['The spectral gap of the graph is ',num2str(spg),'.']);

%% reconstruct signal

tic
xhat = reconstructPolarized(W,Wphase,y,phi,fracRmv,spg);
xhat = relativeOptimalPhase(x,xhat)*xhat; % up to global shift
t = toc;

fprintf('\n');
disp(['Elapsed Time = ', num2str(t)]);
fprintf('\n');
RelativeError = norm(x'-xhat')/norm(x); % L2 error
disp(['Relative Error (L2) = ', num2str(RelativeError)]);
fprintf('\n');

end