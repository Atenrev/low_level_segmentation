function [edgePot,edgeStruct]=CreateGridUGMModel(NumFils, NumCols, K, lambda)
%
%
% NumFils, NumCols: image dimension
% K: number of states
% lambda: smoothing factor

nRows = NumFils;
nCols = NumCols;
nNodes = nRows*nCols;
nStates = K; %2
 
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

% Standardize Features
Xstd = UGM_standardizeCols(reshape(X,[1 1 nNodes]),1);

%The last line transforms the noisy image intensities so that they have a mean of zero and a standard deviation of 1. We will use the following node potentials:
nodePot = zeros(nNodes,nStates);
nodePot(:,1) = exp(-1-2.5*Xstd(:));
nodePot(:,2) = 1;

% We want to use edge potentials that reflect that neighboring pixels are more likely to have the same label. Further, we expect that they are even more likely to have the same label if the difference in their intensities is small. We incoporate this intuition by using the following Ising-like edge potentials
edgePot = zeros(nStates,nStates,edgeStruct.nEdges);

for e = 1:edgeStruct.nEdges
   n1 = edgeStruct.edgeEnds(e,1);
   n2 = edgeStruct.edgeEnds(e,2);

   pot_same = exp(1.8 + .3*1/(1+abs(Xstd(n1)-Xstd(n2))));
   edgePot(:,:,e) = [pot_same 1;1 pot_same];
end
