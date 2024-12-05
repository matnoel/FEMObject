function plotPartition(S,varargin)
% function plotPartition(S)
% Display partition of model S
% S: MODEL

p = ImprovedInputParser;
addParameter(p,'legend',true,@islogical);
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

figure('Name','Mesh partition')
% set(gcf,'Name','Mesh partition')
clf
plotparamelem(S,'partition');
if ~p.Results.legend
    legend('off')
end
set(gca,'FontSize',p.Results.FontSize)

end