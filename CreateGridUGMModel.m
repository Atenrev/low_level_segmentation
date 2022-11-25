function [edgePot,edgeStruct]=CreateGridUGMModel(NumFils, NumCols, K, lambda)
%
%
% NumFils, NumCols: image dimension
% K: number of states
% lambda: smoothing factor

nRows = NumFils;
nCols = NumCols;
nNodes = nRows*nCols;
nStates = K;
 
adj = sparse(nNodes,nNodes);
 
% Add Down Edges
ind = 1:nNodes;
exclude = sub2ind([nRows nCols],repmat(nRows,[1 nCols]),1:nCols); % No Down edge for last row
ind = setdiff(ind,exclude);
adj(sub2ind([nNodes nNodes],ind,ind+1)) = 1;
 
% Add Right Edges
ind = 1:nNodes;
exclude = sub2ind([nRows nCols],1:nRows,repmat(nCols,[1 nRows])); % No right edge for last column
ind = setdiff(ind,exclude);
adj(sub2ind([nNodes nNodes],ind,ind+nRows)) = 1;
 
% Add Up/Left Edges
adj = adj+adj';
edgeStruct = UGM_makeEdgeStruct(adj,nStates);

% We want to use edge potentials that reflect that neighboring pixels are more likely to have the same label. Further, we expect that they are even more likely to have the same label if the difference in their intensities is small. We incoporate this intuition by using the following Ising-like edge potentials
edgePot = zeros(nStates,nStates,edgeStruct.nEdges);

for e = 1:edgeStruct.nEdges
   %n1 = edgeStruct.edgeEnds(e,1);
   %n2 = edgeStruct.edgeEnds(e,2);

   %pot_same = exp(1.8 + .3*1/(1+abs(Xstd(n1)-Xstd(n2))));
   zero_diag = ones(K,K)-eye(K); % 
   %pots = exp(-zero_diag*lambda); % Potencial de parelles exponencial
   %pots = exp(-zero_diag * (rand(K, K)*lambda));
   pots = zero_diag * rand(K, K)*lambda;
   edgePot(:,:,e) = pots; % Inicialitzem edges amb Pots 
end
