function plotAllSolutionsReference(glob,patches,interfaces,U_ref,w_ref,lambda_ref,varargin)
% function plotAllSolutionsReference(glob,patches,interfaces,U_ref,w_ref,lambda_ref,varargin)
% Display reference multiscale solution u_ref=(U_ref,w_ref), global solution U_ref,
% local solution w_ref and Lagrange multiplier lambda_ref
% glob: Global or GlobalOutside
% patches: Patches
% interfaces: Interfaces
% U_ref: reference global solution U
% w_ref: reference local solution w
% lambda_ref: reference Lagrange multiplier lambda

p = ImprovedInputParser;
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

varargin = delcharin({'colorbar','colormap','FontSize'},varargin);
n = numel(patches);

if ischarin('sigma',varargin)
    i = getcharin('sigma',varargin);
    figure('Name',['Reference solution sig_u_ref_' num2str(i) '=(sig_U_ref_' num2str(i) ',sig_w_ref_' num2str(i) '), sig_U_ref_' num2str(i) ', sig_w_ref_' num2str(i) ', sig_lambda_ref_' num2str(i)])
    % set(gcf,'Name',['Reference solution sig_u_ref_' num2str(i) '=(sig_U_ref_' num2str(i) ',sig_w_ref_' num2str(i) '), sig_U_ref_' num2str(i) ', sig_w_ref_' num2str(i) ', sig_lambda_ref_' num2str(i)]])
elseif ischarin('epsilon',varargin)
    i = getcharin('epsilon',varargin);
    figure('Name',['Reference solution eps_u_ref_' num2str(i) '=(eps_U_ref_' num2str(i) ',eps_w_ref_' num2str(i) '), eps_U_ref_' num2str(i) ', eps_w_ref_' num2str(i) ', eps_lambda_ref_' num2str(i)])
    % set(gcf,'Name',['Reference solution eps_u_ref_' num2str(i) '=(eps_U_ref_' num2str(i) ',eps_w_ref_' num2str(i) '), eps_U_ref_' num2str(i) ', eps_w_ref_' num2str(i) ', eps_lambda_ref_' num2str(i)]])
elseif ischarin('energyint',varargin)
    figure('Name',['Reference solution H_u_ref=(H_U_ref,H_w_ref), H_U_ref, H_w_ref, H_lambda_ref'])
    % set(gcf,'Name',['Reference solution H_u_ref=(H_U_ref,H_w_ref), H_U_ref, H_w_ref, H_lambda_ref'])
elseif ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Reference solution u_ref_' num2str(i) '=(U_ref_' num2str(i) ',w_ref_' num2str(i) '), U_ref_' num2str(i) ', w_ref_' num2str(i) ', lambda_ref_' num2str(i)])
    % set(gcf,'Name',['Reference solution u_ref_' num2str(i) '=(U_ref_' num2str(i) ',w_ref_' num2str(i) '), U_ref_' num2str(i) ', w_ref_' num2str(i) ', lambda_ref_' num2str(i)]])
elseif ischarin('rotation',varargin)
    i = getcharin('rotation',varargin);
    figure('Name',['Reference solution rot_u_ref_' num2str(i) '=(rot_U_ref_' num2str(i) ',rot_w_ref_' num2str(i) '), rot_U_ref_' num2str(i) ', rot_w_ref_' num2str(i) ', rot_lambda_ref_' num2str(i)])
    % set(gcf,'Name',['Reference solution rot_u_ref_' num2str(i) '=(rot_U_ref_' num2str(i) ',rot_w_ref_' num2str(i) '), rot_U_ref_' num2str(i) ', rot_w_ref_' num2str(i) ', rot_lambda_ref_' num2str(i)]])
else
    figure('Name','Reference solution u_ref=(U_ref,w_ref), U_ref, w_ref, lambda_ref')
    % set(gcf,'Name','Reference solution u_ref=(U_ref,w_ref), U_ref, w_ref, lambda_ref')
end
clf

if isa(glob,'Global')
    S_out = glob.S_out;
elseif isa(glob,'GlobalOutside')
    S_out = glob.S;
end

subplot(1+n,2,1)
plot_sol(S_out,U_ref,varargin{:});
for k=1:n
    patch = patches.patches{k};
    interface = interfaces.interfaces{k};
    plot_sol(patch.S,w_ref{patch.number},varargin{:});
    if ~ischarin('sigma',varargin) && ~ischarin('epsilon',varargin) && ~ischarin('energyint',varargin)
        plot_sol(interface.S,interface.P_patch*w_ref{patch.number},'FaceColor','none','EdgeColor','k',varargin{:});
    end
end
colormap(p.Results.colormap)
if p.Results.colorbar
    colorbar
end
ax=axis;
cax=caxis;
set(gca,'FontSize',p.Results.FontSize)

subplot(1+n,2,2)
plot_sol(S_out,U_ref,varargin{:});
colormap(p.Results.colormap)
if p.Results.colorbar
    colorbar
end
axis(ax)
% caxis(cax)
set(gca,'FontSize',p.Results.FontSize)

for k=1:n
    subplot(1+n,2,2*k+1)
    patch = patches.patches{k};
    plot_sol(patch.S,w_ref{patch.number},varargin{:});
    colormap(p.Results.colormap)
    if p.Results.colorbar
        colorbar
    end
    % axis(ax)
    % caxis(cax)
    set(gca,'FontSize',p.Results.FontSize)
    
    if ~ischarin('sigma',varargin) && ~ischarin('epsilon',varargin) && ~ischarin('energyint',varargin)
        subplot(1+n,2,2*(k+1))
        interface = interfaces.interfaces{k};
        plot_sol(interface.S,lambda_ref{interface.number},'EdgeColor','interp',varargin{:});
        colormap(p.Results.colormap)
        if p.Results.colorbar
            colorbar
        end
        % axis(ax)
        % caxis(cax)
        set(gca,'FontSize',p.Results.FontSize)
    end
end

end
