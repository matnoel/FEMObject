function plotStd(S,u,varargin)
% function plotStd(S,u,varargin)
% Display the standard deviation of solution u
% S: MODEL
% u: FunctionalBasisArray of solution u

p = ImprovedInputParser;
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

varargin = delcharin({'colorbar','colormap','FontSize'},varargin);

if ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Standard deviation of u_' num2str(i)])
    % set(gcf,'Name',['Standard deviation of u_' num2str(i)])
else
    figure('Name','Standard deviation of u')
    % set(gcf,'Name','Standard deviation of u')
end
clf

s = std(u)';
s = unfreevector(S,s)-calc_init_dirichlet(S);
plot_sol(S,s,varargin{:});
colormap(p.Results.colormap)
if p.Results.colorbar
    colorbar
end
set(gca,'FontSize',p.Results.FontSize)

end
