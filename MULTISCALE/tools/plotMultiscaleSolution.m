function plotMultiscaleSolution(glob,patches,interfaces,U,w,varargin)
% function plotMultiscaleSolution(glob,patches,interfaces,U,w,varargin)
% Display multiscale solution u=(U,w)
% glob: Global
% patches: Patches
% interfaces: Interfaces
% U: global solution U
% w: local solution w

p = ImprovedInputParser;
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

varargin = delcharin({'colorbar','colormap','FontSize'},varargin);
n = numel(patches);

if ischarin('sigma',varargin)
    i = getcharin('sigma',varargin);
    figure('Name',['Multiscale solution sig_u_' num2str(i) '=(sig_U_' num2str(i) ',sig_w_' num2str(i) ') over domain'])
    % set(gcf,'Name',['Multiscale solution sig_u_' num2str(i) '=(sig_U_' num2str(i) ',sig_w_' num2str(i) ') over domain'])
elseif ischarin('epsilon',varargin)
    i = getcharin('epsilon',varargin);
    figure('Name',['Multiscale solution eps_u_' num2str(i) '=(eps_U_' num2str(i) ',eps_w_' num2str(i) ') over domain'])
    % set(gcf,'Name',['Multiscale solution eps_u_' num2str(i) '=(eps_U_' num2str(i) ',eps_w_' num2str(i) ') over domain'])
elseif ischarin('energyint',varargin)
    figure('Name',['Multiscale solution H_u=(H_U,H_w) over domain'])
    % set(gcf,'Name',['Multiscale solution H_u=(H_U,H_w) over domain'])
elseif ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Multiscale solution u_' num2str(i) '=(U_' num2str(i) ',w_' num2str(i) ') over domain'])
    % set(gcf,'Name',['Multiscale solution u_' num2str(i) '=(U_' num2str(i) ',w_' num2str(i) ') over domain'])
elseif ischarin('rotation',varargin)
    i = getcharin('rotation',varargin);
    figure('Name',['Multiscale solution rot_u_' num2str(i) '=(rot_U_' num2str(i) ',rot_w_' num2str(i) ') over domain'])
    % set(gcf,'Name',['Multiscale solution rot_u_' num2str(i) '=(rot_U_' num2str(i) ',rot_w_' num2str(i) ') over domain'])
else
    figure('Name','Multiscale solution u=(U,w) over domain')
    % set(gcf,'Name','Multiscale solution u=(U,w) over domain')
end
clf

if size(U,1)==glob.S.nbddl
    P_out = calcProjection(glob,'free',false);
else
    P_out = glob.P_out;
end

plot_sol(glob.S_out,P_out*U,varargin{:});
for k=1:n
    patch = patches.patches{k};
    interface = interfaces.interfaces{k};
    plot_sol(patch.S,w{patch.number},varargin{:});
    if ~ischarin('sigma',varargin) && ~ischarin('epsilon',varargin) && ~ischarin('energyint',varargin)
        plot_sol(interface.S,interface.P_patch*w{patch.number},'FaceColor','none','EdgeColor','k',varargin{:});
    end
end
colormap(p.Results.colormap)
if p.Results.colorbar
    colorbar
end
set(gca,'FontSize',p.Results.FontSize)

end
