function patches = batch(varargin)
% patches = batch(varargin)
% Create multiple instances of QPPatch at once and return them in a cell
% array. 
% Input arguments are the same as for QPPatch constructor: a list of pairs
% "name" and "value"; the difference is that here "value" is a
% cell array of such values as expected by QPPatch constructor. All those
% arrays must have same length, which will be output argument length.

% Get input arguments
globalModel = getcharin('model',varargin) ;
types = getcharin('type',varargin) ;
K = getcharin('conductivity',varargin) ;
solvers = getcharin('solver',varargin) ;
cells = getcharin('cells',varargin) ;
% Uniformize cells input argument
if iscell(cells{1}) % Keep only first cell list
    cells = cellfun(@(x) x{1},cells(:),'UniformOutput',false) ;
end
patchesNb = numel(types) ;
% In case of unique type or solver
if numel(types) == 1
    types = repmat(types,1,patchesNb) ;
end
if numel(solvers) == 1
    solvers = repmat(solvers,1,patchesNb) ;
end
% Create patches
patches = cell(patchesNb,1) ;
for n = 1:patchesNb
    pModel = subModel(globalModel,cells{n}) ;
    patches{n} = QPPatch('model',pModel,...
        'type',types(n),...
        'conductivity',K{n},...
        'solver',solvers(n)) ;
end
end