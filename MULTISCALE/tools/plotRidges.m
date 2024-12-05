function plotRidges(S,varargin)
% function plotRidges(S)
% Display ridges of model S
% S: MODEL

p = ImprovedInputParser;
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

figure('Name','Mesh ridges')
% set(gcf,'Name','Mesh ridges')
clf
h = plotridges(S);
hg = hggroup;
set(h(:),'Parent',hg);
set(gca,'FontSize',p.Results.FontSize)

end