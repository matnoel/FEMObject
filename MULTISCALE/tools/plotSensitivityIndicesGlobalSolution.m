function plotSensitivityIndicesGlobalSolution(glob,U,alpha,varargin)
% function plotSensitivityIndicesGlobalSolution(glob,U,alpha,varargin)
% Display the sensitivity indices of global solution U associated with the
% group of variables alpha in {1,..,d} with d = ndims(U)
% glob: Global
% U: FunctionalBasisArray of global solution U
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
    figure('Name',['Sensitivity index of U_' num2str(i) ' over fictitious domain for random variables #' num2str(alpha)])
    % set(gcf,'Name',['Sensitivity index of U_' num2str(i) ' over fictitious domain for random variables #' num2str(alpha)])
else
    figure('Name',['Sensitivity index of U over fictitious domain for random variables #' num2str(alpha)])
    % set(gcf,'Name',['Sensitivity index of U over fictitious domain for random variables #' num2str(alpha)])
end
clf

v = variance(U);
s = varianceConditionalExpectation(U,alpha)'./max(v);
s = unfreevector(glob.S,s)-calc_init_dirichlet(glob.S);
plot_sol(glob.S,s,varargin{:});
colormap(p.Results.colormap)
if p.Results.colorbar
    colorbar
end
set(gca,'FontSize',p.Results.FontSize)

end
