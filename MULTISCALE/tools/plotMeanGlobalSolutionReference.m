function plotMeanGlobalSolutionReference(glob,U_ref,varargin)
% function plotMeanGlobalSolutionReference(glob,U_ref,varargin)
% Display the mean (mathematical expectation) of reference global solution U_ref
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
    figure('Name',['Mean of U_ref_' num2str(i) ' over complementary subdomain'])
    % set(gcf,'Name',['Mean of U_ref_' num2str(i) ' over complementary subdomain'])
else
    figure('Name','Mean of U_ref over complementary subdomain')
    % set(gcf,'Name','Mean of U_ref over complementary subdomain')
end
clf

if isa(glob,'Global')
    S_out = glob.S_out;
elseif isa(glob,'GlobalOutside')
    S_out = glob.S;
end

plot_sol(S_out,mean(U_ref)',varargin{:});
colormap(p.Results.colormap)
if p.Results.colorbar
    colorbar
end
set(gca,'FontSize',p.Results.FontSize)

end
