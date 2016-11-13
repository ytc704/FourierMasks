function [W,Wphase,y,phi,spg,newMu] = makeMeasurements(x,k,VarNoise)
% makeMeasurements Returns:
%   'W' adjacency matrix
%   'Wphase' relative phase matrix
%   'phi' original frames
%   'spg' spectral gap
%   'newMu' orignal and additional masks.
%
% Makes noisy measurements with parameter 'VarNoise' using 'k' random masks
% plus additional masks.
%
% Yutong Chen, Princeton University
% Afonso Bandeira, Princeton University
% Dustin Mixon, Air Force Institute of Technology
%
% Free to use. Please cite our paper "Phase retrieval from power spectra of
% masked signals" if the code is used in publications.

M = length(x);

% minimum spectral gap
spgTol = 0.2; % spectral gap tolerance
spg = 0;
while spg < spgTol
%% original masks

mu = random('normal',0,1,M,k);

% frames (goes through all the frequencies first on each mask)
DFT = dftmtx(M)/sqrt(M);
phi = zeros(M,k*M);
for i = 1:M
    phi(:,(1+(i-1)*k):(i*k)) = mu.*repmat(transpose(DFT(i,:)),1,k);
end

% measurements
y = abs(phi'*conj(x)).^2 + random('normal',0,sqrt(VarNoise),k*M,1);

%% additional masks

% Acalc: vector with the indices of the elements in the set \AAA (refer to paper)
sizeAcalc = 0;
while sizeAcalc == 0
    ZM = rand(M-1,1);
    Acalc = [];
    for i = 1:(M-1)
        if ZM(i) <= log(M)/M
            Acalc = [Acalc i];
        end
    end
    Acalc = unique(mod(union(Acalc , -Acalc),M));
    sizeAcalc = length(Acalc);
end

% build adjacency and (not) relative phase matrices - CAN BE MADE MUCH FASTER WITH FFT
w = complex(cos(2*pi/3),sin(2*pi/3));
W = zeros(k*M); % adjencency
Wphase = zeros(k*M); % (not) relative phase --> WILL BE NORMALIZED IN RECONSTRUCTION
newMu = mu;
for r1 = 1:k
for r2 = r1:k
    for i = 1:length(Acalc)
    
        alpha = Acalc(i);
        masks = getMasks(mu,r1,r2,alpha); % computes the 3 polarized vectors
        newMu = [newMu masks'];
    
    for j = 1:M
    
        % edge (m,n)
        m = (j-1)*k + r1; % index of first vector in the polarization
        n = mod(((j-1)-alpha),M)*k + r2; % index of second vector in the polarization
    
        % adjacency
        W(m,n) = 1;
        W(n,m) = W(m,n);
    
        % relative phase
        Wphase(n,m) = sum([w^0 w^1 w^2]*abs((masks.*repmat(DFT(j,:),3,1))*conj(x)).^2)/3 + random('normal',0,sqrt(VarNoise)); % inner products times conjugate of the other (refer to paper)
        Wphase(n,m) = Wphase(n,m);
        Wphase(m,n) = conj(Wphase(n,m));
    
    end
    end
end
end

%% check spectral gap

lambda1 = max(eig(W));
[val,~] = sort(abs(eig(W)),'descend');
lambda2 = val(2);
spg = (lambda1 - lambda2)/(k*length(Acalc));

end

end