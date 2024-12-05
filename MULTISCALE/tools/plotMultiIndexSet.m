function plotMultiIndexSet(f,varargin)
% function plotMultiIndexSet(f,'dim',dim)
% Plot the multi-index set of FunctionalBasisArray f projected on dimension dim
% dim = 1:min(3,PC.M) by default
% 
% function plotMultiIndexSet(f,'dim',dim,'MarkerArea',MarkerArea,'MarkerColor',MarkerColor,'MarkerType',MarkerType,'axis',axis)
% Plot the multi-index set of FunctionalBasisArray f projected on dimension dim with
% optional parameters
% marker = 'o' by default
% color = 'b' by default
% 
% function plotMultiIndexSet(f,'dim',dim,'maximal',true)
% Plot the multi-index set of FunctionalBasisArray f projected on dimension dim
% and its maximal indices
% 
% function plotMultiIndexSet(f,'dim',dim,'margin',true)
% Plot the multi-index set of FunctionalBasisArray f projected on dimension dim
% and its margin
% 
% function plotMultiIndexSet(f,'dim',dim,'reducedmargin',true)
% Plot the multi-index set of FunctionalBasisArray f projected on dimension dim
% and its reduced margin

figure('Name','Multi-index set')
clf
plotMultiIndices(f,varargin{:})

end
