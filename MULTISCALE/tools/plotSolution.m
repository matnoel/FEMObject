function plotSolution(S,u,varargin)
% function plotSolution(S,u,varargin)
% Display solution u
% S: MODEL
% u: solution u

p = ImprovedInputParser;
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

varargin = delcharin({'colorbar','colormap','FontSize'},varargin);

if ischarin('sigma',varargin)
    i = getcharin('sigma',varargin);
    figure('Name',['Solution sig_' num2str(i)])
    % set(gcf,'Name',['Solution sig_' num2str(i)])
elseif ischarin('epsilon',varargin)
    i = getcharin('epsilon',varargin);
    figure('Name',['Solution eps_' num2str(i)])
    % set(gcf,'Name',['Solution eps_' num2str(i)])
elseif ischarin('energyint',varargin)
    figure('Name',['Solution H'])
    % set(gcf,'Name',['Solution H'])
elseif ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Solution u_' num2str(i)])
    % set(gcf,'Name',['Solution u_' num2str(i)])
elseif ischarin('rotation',varargin)
    i = getcharin('rotation',varargin);
    figure('Name',['Solution r_' num2str(i)])
    % set(gcf,'Name',['Solution r_' num2str(i)])
else
    figure('Name','Solution')
    % set(gcf,'Name','Solution')
end
clf

plot_sol(S,u,varargin{:});
colormap(p.Results.colormap)
if p.Results.colorbar
    colorbar
end
set(gca,'FontSize',p.Results.FontSize)

end
