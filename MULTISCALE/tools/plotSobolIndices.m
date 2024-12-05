function plotSobolIndices(S,u,alpha,varargin)
% function plotSobolIndices(S,u,alpha,varargin)
% Display the Sobol indices of solution u associated with the
% group of variables alpha in {1,..,d} with d = ndims(u)
% S: MODEL
% u: FunctionalBasisArray of solution u
% alpha: 1-by-s array of integers or 1-by-d logical
% - if alpha is an array of integers, indices with respect
% to variables alpha
% - if alpha is logical, indices with respect
% to variables find(alpha)

p = ImprovedInputParser;
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

varargin = delcharin({'colorbar','colormap','FontSize'},varargin);

if ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Sobol index of u_' num2str(i) ' for random variables #' num2str(alpha)])
    % set(gcf,'Name',['Sobol index of u_' num2str(i) ' for random variables #' num2str(alpha)])
else
    figure('Name',['Sobol index of u for random variables #' num2str(alpha)])
    % set(gcf,'Name',['Sobol index of u for random variables #' num2str(alpha)])
end
clf

d = ndims(u);
s = SensitivityAnalysis.sobolIndices(u,alpha,d)';
s = unfreevector(S,s)-calc_init_dirichlet(S);
plot_sol(S,s,varargin{:});
colormap(p.Results.colormap)
if p.Results.colorbar
    colorbar
end
set(gca,'FontSize',p.Results.FontSize)

end
