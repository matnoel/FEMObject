function plotStdGlobalSolution(glob,U,varargin)
% function plotStdGlobalSolution(glob,U,varargin)
% Display the standard deviation of global solution U
% glob: Global
% U: FunctionalBasisArray of global solution U

p = ImprovedInputParser;
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

varargin = delcharin({'colorbar','colormap','FontSize'},varargin);

if ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Standard deviation of U_' num2str(i) ' over fictitious domain'])
    % set(gcf,'Name',['Standard deviation of U_' num2str(i) ' over fictitious domain'])
else
    figure('Name','Standard deviation of U over fictitious domain')
    % set(gcf,'Name','Standard deviation of U over fictitious domain')
end
clf

sU = std(U)';
sU = unfreevector(glob.S,sU)-calc_init_dirichlet(glob.S);
plot_sol(glob.S,sU,varargin{:});
colormap(p.Results.colormap)
if p.Results.colorbar
    colorbar
end
set(gca,'FontSize',p.Results.FontSize)

end
