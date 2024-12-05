function plotGlobalSolution(glob,U,varargin)
% function plotGlobalSolution(glob,U,varargin)
% Display global solution U
% glob: Global
% U: global solution U

p = ImprovedInputParser;
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

varargin = delcharin({'colorbar','colormap','FontSize'},varargin);

if ischarin('sigma',varargin)
    i = getcharin('sigma',varargin);
    figure('Name',['Global solution sig_U_' num2str(i) ' over fictitious domain'])
    % set(gcf,'Name',['Global solution sig_U_' num2str(i) ' over fictitious domain'])
elseif ischarin('epsilon',varargin)
    i = getcharin('epsilon',varargin);
    figure('Name',['Global solution eps_U_' num2str(i) ' over fictitious domain'])
    % set(gcf,'Name',['Global solution eps_U_' num2str(i) ' over fictitious domain'])
elseif ischarin('energyint',varargin)
    figure('Name',['Global solution H_U over fictitious domain'])
    % set(gcf,'Name',['Global solution H_U over fictitious domain'])
elseif ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Global solution U_' num2str(i) ' over fictitious domain'])
    % set(gcf,'Name',['Global solution U_' num2str(i) ' over fictitious domain'])
elseif ischarin('rotation',varargin)
    i = getcharin('rotation',varargin);
    figure('Name',['Global solution rot_U_' num2str(i) ' over fictitious domain'])
    % set(gcf,'Name',['Global solution rot_U_' num2str(i) ' over fictitious domain'])
else
    figure('Name','Global solution U over fictitious domain')
    % set(gcf,'Name','Global solution U over fictitious domain')
end
clf

plot_sol(glob.S,U,varargin{:});
colormap(p.Results.colormap)
if p.Results.colorbar
    colorbar
end
set(gca,'FontSize',p.Results.FontSize)

end
