function optimalphase = relativeOptimalPhase(x,u)
% relativeOptimalPhase Returns global shift constant 'optimalphase' between
% 'x' and 'u'.
%
% Yutong Chen, Princeton University
% Afonso Bandeira, Princeton University
% Dustin Mixon, Air Force Institute of Technology
%
% Free to use. Please cite our paper "Phase retrieval from power spectra of
% masked signals" if the code is used in publications.

% w*u = x
w = u'*(conj(x));
w = w/abs(w);

optimalphase = w;

end