function plot_multi_index_set(PC,varargin)
% function plot_multi_index_set(PC,'dim',dim)
% Plot the multi-index set of POLYCHAOS PC projected on dimension dim
% dim = 1:min(3,PC.M) by default
% 
% function plot_multi_index_set(PC,'dim',dim,'MarkerArea',MarkerArea,'MarkerColor',MarkerColor,'MarkerType',MarkerType,'axis',axis)
% Plot the multi-index set of POLYCHAOS PC projected on dimension dim with
% optional parameters
% marker = 'o' by default
% color = 'b' by default
% 
% function plot_multi_index_set(PC,'dim',dim,'maximal',true)
% Plot the multi-index set of POLYCHAOS PC projected on dimension dim
% and its maximal indices
% 
% function plot_multi_index_set(PC,'dim',dim,'margin',true)
% Plot the multi-index set of POLYCHAOS PC projected on dimension dim
% and its margin
% 
% function plot_multi_index_set(PC,'dim',dim,'reducedmargin',true)
% Plot the multi-index set of POLYCHAOS PC projected on dimension dim
% and its reduced margin

figure('Name','Multi-index set')
clf
plot_indices(PC,varargin{:})

end
