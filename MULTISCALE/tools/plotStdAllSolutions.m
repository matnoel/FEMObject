function plotStdAllSolutions(glob,patches,interfaces,U,w,lambda,varargin)
% function plotStdAllSolutions(glob,patches,interfaces,U,w,lambda,varargin)
% Display the standard deviation of multiscale solution u=(U,w),
% global solution U, local solution w and Lagrange multiplier lambda
% glob: Global
% patches: Patches
% interfaces: Interfaces
% U: FunctionalBasisArray of global solution U
% w: FunctionalBasisArray of local solution w
% lambda: FunctionalBasisArray of Lagrange multiplier lambda

p = ImprovedInputParser;
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

varargin = delcharin({'colorbar','colormap','FontSize'},varargin);
n = numel(patches);

if ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Standard deviation of u_' num2str(i) '=(U_' num2str(i) ',w_' num2str(i) '), U_' num2str(i) ', w_' num2str(i) ', lambda_' num2str(i)])
    % set(gcf,'Name',['Standard deviation of u_' num2str(i) '=(U_' num2str(i) ',w_' num2str(i) '), U_' num2str(i) ', w_' num2str(i) ', lambda_' num2str(i)])
else
    figure('Name','Standard deviation of u=(U,w), U, w, lambda')
    % set(gcf,'Name','Standard deviation of u=(U,w), U, w, lambda')
end
clf

if U.sz(1)==glob.S.nbddl
    P_out = calcProjection(glob,'free',false);
else
    P_out = glob.P_out;
end

subplot(1+n,2,1)
sU_out = std(U*P_out')';
sU_out = unfreevector(glob.S_out,sU_out)-calc_init_dirichlet(glob.S_out);
plot_sol(glob.S_out,sU_out,varargin{:});
for k=1:n
    patch = patches.patches{k};
    interface = interfaces.interfaces{k};
    plot_sol(patch.S,std(w{patch.number})',varargin{:});
    if ~ischarin('sigma',varargin) && ~ischarin('epsilon',varargin) && ~ischarin('energyint',varargin)
        plot_sol(interface.S,std(w{patch.number}*interface.P_patch')','FaceColor','none','EdgeColor','k',varargin{:});
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
sU = std(U)';
sU = unfreevector(glob.S,sU)-calc_init_dirichlet(glob.S);
plot_sol(glob.S,sU,varargin{:});
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
    plot_sol(patch.S,std(w{patch.number})',varargin{:});
    colormap(p.Results.colormap)
    if p.Results.colorbar
        colorbar
    end
    % axis(ax)
    % caxis(cax)
    set(gca,'FontSize',p.Results.FontSize)
    
    subplot(1+n,2,2*(k+1))
    interface = interfaces.interfaces{k};
    plot_sol(interface.S,std(lambda{interface.number})','EdgeColor','interp',varargin{:});
    colormap(p.Results.colormap)
    if p.Results.colorbar
        colorbar
    end
    % axis(ax)
    % caxis(cax)
    set(gca,'FontSize',p.Results.FontSize)
end

end
