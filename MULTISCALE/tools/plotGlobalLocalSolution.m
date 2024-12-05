function plotGlobalLocalSolution(glob,patches,interfaces,U,w,varargin)
% function plotGlobalLocalSolution(glob,patches,interfaces,U,w,varargin)
% Display global solution U and local solution w
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
    figure('Name',['Global solution sig_U_' num2str(i) ' and local solution sig_w_' num2str(i) ' over fictitious domain and patches'])
    % set(gcf,'Name',['Global solution sig_U_' num2str(i) ' and local solution sig_w_' num2str(i) ' over fictitious domain and patches'])
elseif ischarin('epsilon',varargin)
    i = getcharin('epsilon',varargin);
    figure('Name',['Global solution eps_U_' num2str(i) ' and local solution eps_w_' num2str(i) ' over fictitious domain and patches'])
    % set(gcf,'Name',['Global solution eps_U_' num2str(i) ' and local solution eps_w_' num2str(i) ' over fictitious domain and patches'])
elseif ischarin('energyint',varargin)
    figure('Name',['Global solution H_U and local solution H_w over fictitious domain and patches'])
    % set(gcf,'Name',['Global solution H_U and local solution H_w over fictitious domain and patches'])
elseif ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Global solution U_' num2str(i) ' and local solution w_' num2str(i) ' over fictitious domain and patches'])
    % set(gcf,'Name',['Global solution U_' num2str(i) ' and local solution w_' num2str(i) ' over fictitious domain and patches'])
elseif ischarin('rotation',varargin)
    i = getcharin('rotation',varargin);
    figure('Name',['Global solution rot_U_' num2str(i) ' and local solution rot_w_' num2str(i) ' over fictitious domain and patches'])
    % set(gcf,'Name',['Global solution rot_U_' num2str(i) ' and local solution rot_w_' num2str(i) ' over fictitious domain and patches'])
else
    figure('Name','Global solution U and local solution w over fictitious domain and patches')
    % set(gcf,'Name','Global solution U and local solution w over fictitious domain and patches')
end
clf

h1 = subplot(2,1,1);
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
    c1 = colorbar;
    posc1 = get(c1,'Position');
    posh1 = get(h1,'Position');
    colorbar('off');
end
ax1 = axis;
cax1 = caxis;
if p.Results.colorbar
    set(h1,'Position',posh1);
end
set(gca,'FontSize',p.Results.FontSize)

h2 = subplot(2,1,2);
plot_sol(glob.S,U,varargin{:});
for k=1:n
    patch = patches.patches{k};
    interface = interfaces.interfaces{k};
    if ~ischarin('sigma',varargin) && ~ischarin('epsilon',varargin) && ~ischarin('energyint',varargin)
        plot_sol(interface.S,interface.P_patch*w{patch.number},'FaceColor','none','EdgeColor','k',varargin{:});
    end
end
colormap(p.Results.colormap)
if p.Results.colorbar
    c2 = colorbar;
    posc2 = get(c2,'Position');
    posh2 = get(h2,'Position');
    posc2(4) = posc1(2)+posc1(4)-posc2(2);
    set(c2,'Position',posc2);
end
ax2 = axis;
cax2 = caxis;
if p.Results.colorbar
    set(h2,'Position',posh2);
end
set(gca,'FontSize',p.Results.FontSize)

ax = ax1;
ax(1:2:end) = min(ax1(1:2:end),ax2(1:2:end));
ax(2:2:end) = max(ax1(2:2:end),ax2(2:2:end));
cax(1) = min(cax1(1),cax2(1));
cax(2) = max(cax1(2),cax2(2));

axis([h1 h2],ax)
caxis(h1,cax)
caxis(h2,cax)

pos1 = get(h1,'Position');
pos2 = get(h2,'Position');
trans = pos1(2)-(pos2(2)+pos2(4));
pos1(2) = pos1(2)-trans;
set(h1,'Position',pos1);

if p.Results.colorbar
    posc2(4) = posc2(4)-trans;
    set(c2,'Position',posc2);
end

end
