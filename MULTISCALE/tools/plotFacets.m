function plotFacets(S,varargin)
% function plotFactes(S)
% Display facets of model S
% S: MODEL

p = ImprovedInputParser;
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

figure('Name','Mesh facets')
% set(gcf,'Name','Mesh facets')
clf
h = plotfacets(S);
hg = hggroup;
set(h(:),'Parent',hg);
set(gca,'FontSize',p.Results.FontSize)

end