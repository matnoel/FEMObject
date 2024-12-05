function plotGlobalSolutionReference(glob,U_ref,varargin)
% function plotGlobalSolutionReference(glob,U_ref,varargin)
% Display reference global solution U_ref
% glob: Global or GlobalOutside
% U_ref: reference global solution U

p = ImprovedInputParser;
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

varargin = delcharin({'colorbar','colormap','FontSize'},varargin);

if ischarin('sigma',varargin)
    i = getcharin('sigma',varargin);
    figure('Name',['Reference global solution sig_U_ref_' num2str(i) ' over complementary subdomain'])
    % set(gcf,'Name',['Reference global solution sig_U_ref_' num2str(i) ' over complementary subdomain'])
elseif ischarin('epsilon',varargin)
    i = getcharin('epsilon',varargin);
    figure('Name',['Reference global solution eps_U_ref_' num2str(i) ' over complementary subdomain'])
    % set(gcf,'Name',['Reference global solution eps_U_ref_' num2str(i) ' over complementary subdomain'])
elseif ischarin('energyint',varargin)
    figure('Name',['Reference global solution H_U_ref over complementary subdomain'])
    % set(gcf,'Name',['Reference global solution H_U_ref over complementary subdomain'])
elseif ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Reference global solution U_ref_' num2str(i) ' over complementary subdomain'])
    % set(gcf,'Name',['Reference global solution U_ref_' num2str(i) ' over complementary subdomain'])
elseif ischarin('rotation',varargin)
    i = getcharin('rotation',varargin);
    figure('Name',['Reference global solution rot_U_ref_' num2str(i) ' over complementary subdomain'])
    % set(gcf,'Name',['Reference global solution rot_U_ref_' num2str(i) ' over complementary subdomain'])
else
    figure('Name','Reference global solution U_ref over complementary subdomain')
    % set(gcf,'Name','Reference global solution U_ref over complementary subdomain')
end
clf

if isa(glob,'Global')
    S_out = glob.S_out;
elseif isa(glob,'GlobalOutside')
    S_out = glob.S;
end

plot_sol(S_out,U_ref,varargin{:});
colormap(p.Results.colormap)
if p.Results.colorbar
    colorbar
end
set(gca,'FontSize',p.Results.FontSize)

end
