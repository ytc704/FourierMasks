function xhat = reconstructPolarized(W,Wphase,y,phi,fracRmv,spg)
% reconstructPolarized Returns reconstructed signal 'xhat' via the
% polarization algorithm.
%
% The following are generated in 'makeMeasurements':
%   'W' adjacency matrix
%   'Wphase' relative phase matrix
%   'y' observed signal
%   'phi' original frames
%
% 'fracRmv' is a parameter specified by the user (in 'demo').
%
% Yutong Chen, Princeton University
% Afonso Bandeira, Princeton University
% Dustin Mixon, Air Force Institute of Technology
%
% Free to use. Please cite our paper "Phase retrieval from power spectra of
% masked signals" if the code is used in publications.

numVertices = length(W); % number of vertices in graph

% normalize entries in 'Wphase' to get actual relative phases
Wphase_normalized = zeros(numVertices);
edges = [];
for i = 1:(numVertices-1)
    for j = (i+1):numVertices
        if W(i,j) ~= 0
           edges = [edges; i j abs(Wphase(i,j))];
           Wphase_normalized(i,j) = Wphase(i,j)/abs(Wphase(i,j));
           Wphase_normalized(j,i) = Wphase(j,i)/abs(Wphase(j,i));
        end
    end
end

%% prun for reliability

edges = sortrows(edges,3);

numRmvdVertices = ceil(fracRmv*numVertices);

% find bad edges
vtxCount = numVertices;
edgeCount = 0;
while (numVertices-numRmvdVertices) < vtxCount
    
    edgeCount = edgeCount + 1;
    rmvdEdges = edges(1:edgeCount,:);

    Wtmp = W;
    for ii = 1:edgeCount
        
        i = rmvdEdges(ii,1);
        j = rmvdEdges(ii,2);
        
        Wtmp(i,j) = 0;
        Wtmp(j,i) = 0;
        
    end
    
    [~,VVV] = largest_component(Wtmp);
    vtxCount = length(find(VVV));
    
end

% remove bad edges/vertices
rmvdEdges = edges(1:edgeCount,:);
for ii = 1:edgeCount
    
    i = rmvdEdges(ii,1);
    j = rmvdEdges(ii,2);
    
    W(i,j) = 0;
    W(j,i) = 0;
    
    Wphase_normalized(i,j) = 0;
    Wphase_normalized(j,i) = 0;
    
end

% update
[WWW,VVV] = largest_component(W);
idx = find(VVV);
DDD = diag(sum(WWW));
WWWphase_normalized = Wphase_normalized(idx,idx);

% check spectral gap
WWW_normalized = DDD^(-1/2)*WWW*DDD^(-1/2);
lambda1 = max(eig(WWW_normalized));
[val,~] = sort(abs(eig(WWW_normalized)),'descend');
lambda2 = val(2);

newSpg = lambda1 - lambda2;

%% prun for connectivity
p = 1 - fracRmv;
q = 1 - 4*fracRmv;
minSpg = (1/8)*(spg - (1 - 2*(q*(1-q)-(1-p))))^2;
while newSpg < minSpg
   
    L = eye(length(WWW)) - DDD^(-1/2)*WWW*DDD^(-1/2);
    [u,~] = eigs(L,2,'SM');
    u = DDD^(-1/2)*u(:,1);
    
    [~,uID] = sort(abs(u),'ascend');
    h = zeros((length(u)-1),1);
    
    % apply pruning algorithm
    for i = 1:length(h)
        S = uID(1:i);
        Scomp = setdiff(uID,S);

        E = 0;
        degS = 0;
        degScomp = 0;
        for ii = 1:length(S)
            degS = degS + sum(WWW(S(ii),:));
            for jj = 1:length(Scomp)
                E = E + WWW(ii,jj);
                degScomp = degScomp + sum(WWW(Scomp(jj),:));
            end
        end
        degScomp = degScomp/length(S);
        
        h(i) = E/min(degS,degScomp);

    end

    [~,hID] = sort(h,'ascend');
    
    % update
    idx = setdiff(uID,uID(1:hID(1)));
    WWW = WWW(idx,idx);
    [WWW,VVV] = largest_component(WWW);
    idx = find(VVV);
    DDD = diag(sum(WWW));
    WWWphase_normalized = Wphase_normalized(idx,idx);

    % check spectral gap
    WWW_normalized = DDD^(-1/2)*WWW*DDD^(-1/2);
    lambda1 = max(eig(WWW_normalized));
    [val,~] = sort(abs(eig(WWW_normalized)),'descend');
    lambda2 = val(2);

    newSpg = lambda1 - lambda2;
    
end

fprintf('\n');
disp(['Before "cleaning", there were ',int2str(numVertices),' vertices.'])
disp(['After "cleaning", there are ',int2str(length(WWW)),' vertices, and the spectral gap is ',num2str(newSpg),'.']);
fprintf('\n');

%% angular synchronization

% synchronize
L = DDD - WWWphase_normalized;

[z,~] = eigs(DDD^(-1/2)*L*DDD^(-1/2),1,'SM');
z = DDD^(-1/2)*z;
z = z./abs(z);

% least squares estimate
xhat = phi(:,idx)'\(abs(sqrt(y(idx))).*z);
xhat = conj(xhat);

end