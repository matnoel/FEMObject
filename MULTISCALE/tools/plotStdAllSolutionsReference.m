function plotStdAllSolutionsReference(glob,patches,interfaces,U_ref,w_ref,lambda_ref,varargin)
% function plotStdAllSolutionsReference(glob,patches,interfaces,U_ref,w_ref,lambda_ref,varargin)
% Display the standard deviation of reference multiscale solution u_ref=(U_ref,w_ref),
% global solution U_ref, local solution w_ref and Lagrange multiplier lambda_ref
% glob: Global or GlobalOutside
% patches: Patches
% interfaces: Interfaces
% U_ref: FunctionalBasisArray of reference global solution U
% w_ref: FunctionalBasisArray of reference local solution w
% lambda_ref: FunctionalBasisArray of reference Lagrange multiplier lambda

p = ImprovedInputParser;
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

varargin = delcharin({'colorbar','colormap','FontSize'},varargin);
n = numel(patches);

if ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Standard deviation of u_ref_' num2str(i) '=(U_ref_' num2str(i) ',w_ref_' num2str(i) '), U_ref_' num2str(i) ', w_ref_' num2str(i) ', lambda_ref_' num2str(i)])
    % set(gcf,'Name',['Standard deviation of u_ref_' num2str(i) '=(U_ref_' num2str(i) ',w_ref_' num2str(i) '), U_ref_' num2str(i) ', w_ref_' num2str(i) ', lambda_ref_' num2str(i)]])
else
    figure('Name','Standard deviation of u_ref=(U_ref,w_ref), U_ref, w_ref, lambda_ref')
    % set(gcf,'Name','Standard deviation of u_ref=(U_ref,w_ref), U_ref, w_ref, lambda_ref')
end
clf

if isa(glob,'Global')
    S_out = glob.S_out;
elseif isa(glob,'GlobalOutside')
    S_out = glob.S;
end

subplot(1+n,2,1)
sU_ref = std(U_ref)';
sU_ref = unfreevector(S_out,sU_ref)-calc_init_dirichlet(S_out);
plot_sol(S_out,sU_ref,varargin{:});
for k=1:n
    patch = patches.patches{k};
    interface = interfaces.interfaces{k};
    plot_sol(patch.S,std(w_ref{patch.number})',varargin{:});
    if ~ischarin('sigma',varargin) && ~ischarin('epsilon',varargin) && ~ischarin('energyint',varargin)
        plot_sol(interface.S,std(w_ref{patch.number}*interface.P_patch')','FaceColor','none','EdgeColor','k',varargin{:});
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
plot_sol(S_out,sU_ref,varargin{:});
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
    plot_sol(patch.S,std(w_ref{patch.number})',varargin{:});
    colormap(p.Results.colormap)
    if p.Results.colorbar
        colorbar
    end
    % axis(ax)
    % caxis(cax)
    set(gca,'FontSize',p.Results.FontSize)
    
    subplot(1+n,2,2*(k+1))
    interface = interfaces.interfaces{k};
    plot_sol(interface.S,std(lambda_ref{interface.number})','EdgeColor','interp',varargin{:});
    colormap(p.Results.colormap)
    if p.Results.colorbar
        colorbar
    end
    % axis(ax)
    % caxis(cax)
    set(gca,'FontSize',p.Results.FontSize)
end

end
