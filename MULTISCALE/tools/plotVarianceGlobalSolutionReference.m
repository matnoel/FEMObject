function plotVarianceGlobalSolutionReference(glob,U_ref,varargin)
% function plotVarianceGlobalSolutionReference(glob,U_ref,varargin)
% Display the variance of reference global solution U_ref
% glob: Global or GlobalOutside
% U_ref: FunctionalBasisArray of reference global solution U

p = ImprovedInputParser;
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

varargin = delcharin({'colorbar','colormap','FontSize'},varargin);

if ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Variance of U_ref_' num2str(i) ' over complementary subdomain'])
    % set(gcf,'Name',['Variance of U_ref_' num2str(i) ' over complementary subdomain'])
else
    figure('Name','Variance of U_ref over complementary subdomain')
    % set(gcf,'Name','Variance of U_ref over complementary subdomain')
end
clf

if isa(glob,'Global')
    S_out = glob.S_out;
elseif isa(glob,'GlobalOutside')
    S_out = glob.S;
end

vU_ref = variance(U_ref)';
vU_ref = unfreevector(S_out,vU_ref)-calc_init_dirichlet(S_out);
plot_sol(S_out,vU_ref,varargin{:});
colormap(p.Results.colormap)
if p.Results.colorbar
    colorbar
end
set(gca,'FontSize',p.Results.FontSize)

end
