function plotVarianceMultiscaleSolutionReference(glob,patches,interfaces,U_ref,w_ref,varargin)
% function plotVarianceMultiscaleSolutionReference(glob,patches,interfaces,U_ref,w_ref,varargin)
% Display the variance of reference multiscale solution u_ref=(U_ref,w_ref)
% glob: Global or GlobalOutside
% patches: Patches
% interfaces: Interfaces
% U_ref: FunctionalBasisArray of reference global solution U
% w_ref: FunctionalBasisArray of reference local solution w

p = ImprovedInputParser;
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

varargin = delcharin({'colorbar','colormap','FontSize'},varargin);
n = numel(patches);

if ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Variance of u_ref_' num2str(i) '=(U_ref_' num2str(i) ',w_ref_' num2str(i) ') over domain'])
    % set(gcf,'Name',['Variance of u_ref_' num2str(i) '=(U_ref_' num2str(i) ',w_ref_' num2str(i) ') over domain'])
else
    figure('Name','Variance of u_ref=(U_ref,w_ref) over domain')
    % set(gcf,'Name','Variance of u_ref=(U_ref,w_ref) over domain')
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
for k=1:n
    patch = patches.patches{k};
    interface = interfaces.interfaces{k};
    plot_sol(patch.S,variance(w_ref{patch.number})',varargin{:});
    if ~ischarin('sigma',varargin) && ~ischarin('epsilon',varargin) && ~ischarin('energyint',varargin)
        plot_sol(interface.S,variance(w_ref{patch.number}*interface.P_patch')','FaceColor','none','EdgeColor','k',varargin{:});
    end
end
colormap(p.Results.colormap)
if p.Results.colorbar
    colorbar
end
set(gca,'FontSize',p.Results.FontSize)

end
