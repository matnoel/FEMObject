function plotVariance(S,u,varargin)
% function plotVariance(S,u,varargin)
% Display the variance of solution u
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
    figure('Name',['Variance of u_' num2str(i)])
    % set(gcf,'Name',['Variance of u_' num2str(i)])
else
    figure('Name','Variance of u')
    % set(gcf,'Name','Variance of u')
end
clf

v = variance(u)';
v = unfreevector(S,v)-calc_init_dirichlet(S);
plot_sol(S,v,varargin{:});
colormap(p.Results.colormap)
if p.Results.colorbar
    colorbar
end
set(gca,'FontSize',p.Results.FontSize)

end
