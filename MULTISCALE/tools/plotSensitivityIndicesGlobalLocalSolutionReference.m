function plotSensitivityIndicesGlobalLocalSolutionReference(glob,patches,interfaces,U_ref,w_ref,alpha,varargin)
% function plotSensitivityIndicesGlobalLocalSolutionReference(glob,patches,interfaces,U_ref,w_ref,alpha,varargin)
% Display the sensitivity indices of reference global solution U_ref and local solution w_ref
% associated with the group of variables alpha in {1,..,d} with d = min(ndims(U_ref),ndims(w_ref))
% glob: Global or GlobalOutside
% patches: Patches
% interfaces: Interfaces
% U_ref: FunctionalBasisArray of reference global solution U
% w_ref: FunctionalBasisArray of reference local solution w
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
    figure('Name',['Sensitivity index of U_ref_' num2str(i) ' and w_ref_' num2str(i) ' over complementary subdomain and patches for random variables #' num2str(alpha)])
    % set(gcf,'Name',['Sensitivity index of U_ref_' num2str(i) ' and w_ref_' num2str(i) ' over complementary subdomain and patches for random variables #' num2str(alpha)])
else
    figure('Name',['Sensitivity index of U_ref and w_ref over complementary subdomain and patches for random variables #' num2str(alpha)])
    % set(gcf,'Name',['Sensitivity index of U_ref and w_ref over complementary subdomain and patches for random variables #' num2str(alpha)])
end
clf

if isa(glob,'Global')
    S_out = glob.S_out;
elseif isa(glob,'GlobalOutside')
    S_out = glob.S;
end

max_var = max(variance(U_ref));
for k=1:n
    patch = patches.patches{k};
    max_var = max(max_var,max(variance(w_ref{patch.number})));
end

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
    sw_ref = varianceConditionalExpectation(w_ref{patch.number},alpha)'./max_var;
    plot_sol(patch.S,sw_ref,varargin{:});
    if ~ischarin('sigma',varargin) && ~ischarin('epsilon',varargin) && ~ischarin('energyint',varargin)
        plot_sol(interface.S,varianceConditionalExpectation(w_ref{patch.number}*interface.P_patch',alpha)'./max_var,'FaceColor','none','EdgeColor','k',varargin{:});
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
sU_ref = varianceConditionalExpectation(U_ref,alpha)'./max_var;
sU_ref = unfreevector(S_out,sU_ref)-calc_init_dirichlet(S_out);
plot_sol(S_out,sU_ref,varargin{:});
for k=1:n
    patch = patches.patches{k};
    interface = interfaces.interfaces{k};
    if ~ischarin('sigma',varargin) && ~ischarin('epsilon',varargin) && ~ischarin('energyint',varargin)
        plot_sol(interface.S,varianceConditionalExpectation(w_ref{patch.number}*interface.P_patch',alpha)'./max_var,'FaceColor','none','EdgeColor','k',varargin{:});
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
