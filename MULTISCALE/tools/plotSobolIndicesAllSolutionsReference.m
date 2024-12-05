function plotSobolIndicesAllSolutionsReference(glob,patches,interfaces,U_ref,w_ref,lambda_ref,alpha,varargin)
% function plotSobolIndicesAllSolutionsReference(glob,patches,interfaces,U_ref,w_ref,lambda_ref,alpha,varargin)
% Display Sobol indices of reference multiscale solution u_ref=(U_ref,w_ref),
% global solution U_ref, local solution w_ref and Lagrange multiplier lambda_ref
% associated with the group of variables alpha in {1,..,d} with d = ndims(u_ref)
% glob: Global
% patches: Patches
% interfaces: Interfaces
% U_ref: FunctionalBasisArray of reference global solution U
% w_ref: FunctionalBasisArray of reference local solution w
% lambda_ref: FunctionalBasisArray of reference Lagrange multiplier lambda
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
n = numel(patches);

if ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Sobol index of u_ref_' num2str(i) '=(U_ref_' num2str(i) ',w_ref_' num2str(i) '), U_ref_' num2str(i) ', w_ref_' num2str(i) ', lambda_ref_' num2str(i) ' for random variables #' num2str(alpha)])
    % set(gcf,'Name',['Sobol index of u_ref_' num2str(i) '=(U_ref_' num2str(i) ',w_ref_' num2str(i) '), U_ref_' num2str(i) ', w_ref_' num2str(i) ', lambda_ref_' num2str(i) ' for random variables #' num2str(alpha)])
else
    figure('Name',['Sobol index of u_ref=(U_ref,w_ref), U_ref, w_ref, lambda_ref for random variables #' num2str(alpha)'])
    % set(gcf,'Name',['Sobol index of u_ref=(U_ref,w_ref), U_ref, w_ref, lambda_ref for random variables #' num2str(alpha)'])
end
clf

if isa(glob,'Global')
    S_out = glob.S_out;
elseif isa(glob,'GlobalOutside')
    S_out = glob.S;
end

subplot(1+n,2,1)
dU_ref = ndims(U_ref);
sU_ref = SensitivityAnalysis.sobolIndices(U_ref,alpha,dU_ref)';
sU_ref = unfreevector(S_out,sU_ref)-calc_init_dirichlet(S_out);
plot_sol(S_out,sU_ref,varargin{:});
for k=1:n
    patch = patches.patches{k};
    interface = interfaces.interfaces{k};
    dw_ref = ndims(w_ref{patch.number});
    plot_sol(patch.S,SensitivityAnalysis.sobolIndices(w_ref{patch.number},alpha,dw_ref)',varargin{:});
    if ~ischarin('sigma',varargin) && ~ischarin('epsilon',varargin) && ~ischarin('energyint',varargin)
        plot_sol(interface.S,SensitivityAnalysis.sobolIndices(w_ref{patch.number}*interface.P_patch',alpha,dw_ref)','FaceColor','none','EdgeColor','k',varargin{:});
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
    dw_ref = ndims(w_ref{patch.number});
    plot_sol(patch.S,SensitivityAnalysis.sobolIndices(w_ref{patch.number},alpha,dw_ref)',varargin{:});
    colormap(p.Results.colormap)
    if p.Results.colorbar
        colorbar
    end
    % axis(ax)
    % caxis(cax)
    set(gca,'FontSize',p.Results.FontSize)
    
    subplot(1+n,2,2*(k+1))
    interface = interfaces.interfaces{k};
    dlambda_ref = ndims(lambda_ref{patch.number});
    plot_sol(interface.S,SensitivityAnalysis.sobolIndices(lambda_ref{interface.number},alpha,dlambda_ref)','EdgeColor','interp',varargin{:});
    colormap(p.Results.colormap)
    if p.Results.colorbar
        colorbar
    end
    % axis(ax)
    % caxis(cax)
    set(gca,'FontSize',p.Results.FontSize)
end

end
