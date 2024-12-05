function plotVarianceMultiscaleSolution(glob,patches,interfaces,U,w,varargin)
% function plotVarianceMultiscaleSolution(glob,patches,interfaces,U,w,varargin)
% Display the variance of multiscale solution u=(U,w)
% glob: Global
% patches: Patches
% interfaces: Interfaces
% U: FunctionalBasisArray of global solution U
% w: FunctionalBasisArray of local solution w

p = ImprovedInputParser;
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

varargin = delcharin({'colorbar','colormap','FontSize'},varargin);
n = numel(patches);

if ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Variance of u_' num2str(i) '=(U_' num2str(i) ',w_' num2str(i) ') over domain'])
    % set(gcf,'Name',['Variance of u_' num2str(i) '=(U_' num2str(i) ',w_' num2str(i) ') over domain'])
else
    figure('Name','Variance of u=(U,w) over domain')
    % set(gcf,'Name','Variance of u=(U,w) over domain')
end
clf

if U.sz(1)==glob.S.nbddl
    P_out = calcProjection(glob,'free',false);
else
    P_out = glob.P_out;
end

vU_out = variance(U*P_out')';
vU_out = unfreevector(glob.S_out,vU_out)-calc_init_dirichlet(glob.S_out);
plot_sol(glob.S_out,vU_out,varargin{:});
for k=1:n
    patch = patches.patches{k};
    interface = interfaces.interfaces{k};
    plot_sol(patch.S,variance(w{patch.number})',varargin{:});
    if ~ischarin('sigma',varargin) && ~ischarin('epsilon',varargin) && ~ischarin('energyint',varargin)
        plot_sol(interface.S,variance(w{patch.number}*interface.P_patch')','FaceColor','none','EdgeColor','k',varargin{:});
    end
end
colormap(p.Results.colormap)
if p.Results.colorbar
    colorbar
end
set(gca,'FontSize',p.Results.FontSize)

end
