function plotTotalSobolIndicesGlobalLocalSolution(glob,patches,interfaces,U,w,alpha,varargin)
% function plotTotalSobolIndicesGlobalLocalSolution(glob,patches,interfaces,U,w,alpha,varargin)
% Display the Total Sobol indices of global solution U and local solution w
% associated with the group of variables alpha in {1,..,d} with d = min(ndims(U),ndims(w))
% glob: Global
% patches: Patches
% interfaces: Interfaces
% U: FunctionalBasisArray of global solution U
% w: FunctionalBasisArray of local solution w
% alpha: 1-by-s array of integers or 1-by-d logical
% - if alpha is an array of integers, indices with respect
% to variables alpha
% - if alpha is logical, indices with respect
% to variables find(alpha)

p = ImprovedInputParser;
addParameter(p,'orientation','vertical',@ischar);
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

varargin = delcharin({'colorbar','colormap','FontSize'},varargin);
n = numel(patches);

if ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Total Sobol index of U_' num2str(i) ' and w_' num2str(i) ' over fictitious domain and patches for random variables #' num2str(alpha)])
    % set(gcf,'Name',['Total Sobol index of U_' num2str(i) ' and w_' num2str(i) ' over fictitious domain and patches for random variables #' num2str(alpha)])
else
    figure('Name',['Total Sobol index of U and w over fictitious domain and patches for random variables #' num2str(alpha)])
    % set(gcf,'Name',['Total Sobol index of U and w over fictitious domain and patches for random variables #' num2str(alpha)])
end
clf

if strcmp(p.Results.orientation,'vertical') || strcmp(p.Results.orientation,'v')
    h1 = subplot(2,1,1);
elseif strcmp(p.Results.orientation,'horizontal') || strcmp(p.Results.orientation,'h')
    h1 = subplot(1,2,2);
else
    error('Wrong parameter ''orientation''')
end
for k=1:n
    patch = patches.patches{k};
    interface = interfaces.interfaces{k};
    dw = ndims(w{patch.number});
    plot_sol(patch.S,SensitivityAnalysis.totalSobolIndices(w{patch.number},alpha,dw)',varargin{:});
    if ~ischarin('sigma',varargin) && ~ischarin('epsilon',varargin) && ~ischarin('energyint',varargin)
        plot_sol(interface.S,SensitivityAnalysis.totalSobolIndices(w{patch.number}*interface.P_patch',alpha,dw)','FaceColor','none','EdgeColor','k',varargin{:});
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

if strcmp(p.Results.orientation,'vertical') || strcmp(p.Results.orientation,'v')
    h2 = subplot(2,1,2);
elseif strcmp(p.Results.orientation,'horizontal') || strcmp(p.Results.orientation,'h')
    h2 = subplot(1,2,1);
else
    error('Wrong parameter ''orientation''')
end
dU = ndims(U);
sU = SensitivityAnalysis.totalSobolIndices(U,alpha,dU)';
sU = unfreevector(glob.S,sU)-calc_init_dirichlet(glob.S);
plot_sol(glob.S,sU,varargin{:});
for k=1:n
    patch = patches.patches{k};
    interface = interfaces.interfaces{k};
    if ~ischarin('sigma',varargin) && ~ischarin('epsilon',varargin) && ~ischarin('energyint',varargin)
        dw = ndims(w{patch.number});
        plot_sol(interface.S,SensitivityAnalysis.totalSobolIndices(w{patch.number}*interface.P_patch',alpha,dw)','FaceColor','none','EdgeColor','k',varargin{:});
    end
end
colormap(p.Results.colormap)
if p.Results.colorbar
    c2 = colorbar;
    posc2 = get(c2,'Position');
    posh2 = get(h2,'Position');
    if strcmp(p.Results.orientation,'vertical') || strcmp(p.Results.orientation,'v')
        posc2(4) = posc1(2)+posc1(4)-posc2(2);
    elseif strcmp(p.Results.orientation,'horizontal') || strcmp(p.Results.orientation,'h')
        posc2(1) = posc1(1);
    end
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
if strcmp(p.Results.orientation,'vertical') || strcmp(p.Results.orientation,'v')
    trans = pos1(2)-(pos2(2)+pos2(4));
    pos1(2) = pos1(2)-trans;
    set(h1,'Position',pos1);
    if p.Results.colorbar
        posc2(4) = posc2(4)-trans;
        set(c2,'Position',posc2);
    end
elseif strcmp(p.Results.orientation,'horizontal') || strcmp(p.Results.orientation,'h')
    trans = pos1(1)-(pos2(1)+pos2(3));
    pos1(1) = pos1(1)-trans;
    set(h1,'Position',pos1);
    if p.Results.colorbar
        posc2(1) = posc2(1)-trans;
        set(c2,'Position',posc2);
    end
end

end
