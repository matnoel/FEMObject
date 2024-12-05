function [] = plotModes(model,tensor,gridSize)
% [] = plotModes(model,tensor,gridSize)

if nargin < 3
    gridSize = [] ;
end

modeNb = tensor.space.dim(end) ;
cellModel = getCellModel(model) ;

if isempty(gridSize)
    gridSize = subplot_format(modeNb) ;
end

for i = 1:modeNb
    subplot(gridSize(1),gridSize(2),i)
    plot_sol(cellModel,tensor.space.spaces{end}(:,i))
end

end