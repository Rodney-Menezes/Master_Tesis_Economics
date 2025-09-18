function [vIndicesAbove,vIndicesBelow,vWeightAbove,vWeightBelow] = compute_weighting_matrices(vMapping,vGrid)

% Get the sizes of the grid and the mapping
[nGrid,~]		= size(vGrid);
[nMapping,~]	= size(vMapping);

% Compute the nearest neighbor of implied point in vGrid
vIndices		= knnsearch(vGrid,vMapping);
vGridIndices	= vGrid(vIndices);

% Find indices and gridpoints above and below nearest neighbor
vIndicesBelow = vIndices; vIndicesAbove = vIndices;
vIndicesBelow(vGridIndices >= vMapping) = vIndicesBelow(vGridIndices >= vMapping) - 1;
    vIndicesBelow(vIndicesBelow <= 1) = 1;          % ensure never go below lower bound
    vGridBelow = vGrid(vIndicesBelow);
vIndicesAbove(vGridIndices <= vMapping) = vIndicesAbove(vGridIndices <= vMapping) + 1;
    vIndicesAbove(vIndicesAbove >= nGrid) = nGrid;	% ensure never go above upper bound
    vGridAbove = vGrid(vIndicesAbove);

% Compute weighting matrices
vWeightBelow = (vGridAbove - vMapping) ./ (vGridAbove - vGridBelow);      % mass of points to put onto point below
    vWeightBelow(vMapping <= vGridBelow) = 1;
    vWeightBelow(vMapping >= vGridAbove) = 0;
vWeightAbove = (vMapping - vGridBelow) ./ (vGridAbove - vGridBelow);      % mass of points to put onto point above
    vWeightAbove(vMapping <= vGridBelow) = 0;
    vWeightAbove(vMapping >= vGridAbove) = 1;
