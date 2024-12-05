function plotSensitivityIndicesGlobalSolutionReference(glob,U_ref,alpha,varargin)
% function plotSensitivityIndicesGlobalSolutionReference(glob,U_ref,alpha,varargin)
% Display the sensitivity indices of reference global solution U_ref associated 
% with the group of variables alpha in {1,..,d} with d = ndims(U_ref)
% glob: Global or GlobalOutside
% U_ref: FunctionalBasisArray of reference global solution U
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
    figure('Name',['Sensitivity index of U_ref_' num2str(i) ' over complementary subdomain for random variables #' num2str(alpha)])
    % set(gcf,'Name',['Sensitivity index of U_ref_' num2str(i) ' over complementary subdomain for random variables #' num2str(alpha)])
else
    figure('Name',['Sensitivity index of U_ref over complementary subdomain for random variables #' num2str(alpha)])
    % set(gcf,'Name',['Sensitivity index of U_ref over complementary subdomain for random variables #' num2str(alpha)])
end
clf

if isa(glob,'Global')
    S_out = glob.S_out;
elseif isa(glob,'GlobalOutside')
    S_out = glob.S;
end

v = variance(U);
s = varianceConditionalExpectation(U,alpha)'./max(v);
s = unfreevector(S_out,s)-calc_init_dirichlet(S_out);
plot_sol(S_out,s,varargin{:});
colormap(p.Results.colormap)
if p.Results.colorbar
    colorbar
end
set(gca,'FontSize',p.Results.FontSize)

end
