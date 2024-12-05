function plotLocalSolution(patches,w,varargin)
% function plotLocalSolution(patches,w,varargin)
% Display local solution w
% patches: Patches or Patch
% w: local solution w

p = ImprovedInputParser;
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

varargin = delcharin({'colorbar','colormap','FontSize'},varargin);

if isa(patches,'Patches')
    numbers = getnumber(patches);
    if ischarin('sigma',varargin)
        i = getcharin('sigma',varargin);
        figure('Name',['Local solution sig_w_' num2str(i) ' over patches #' num2str([numbers{:}])])
        % set(gcf,'Name',['Local solution sig_w_' num2str(i) ' over patches #' num2str([numbers{:}])])
    elseif ischarin('epsilon',varargin)
        i = getcharin('epsilon',varargin);
        figure('Name',['Local solution eps_w_' num2str(i) ' over patches #' num2str([numbers{:}])])
        % set(gcf,'Name',['Local solution eps_w_' num2str(i) ' over patches #' num2str([numbers{:}])])
    elseif ischarin('energyint',varargin)
        figure('Name',['Local solution H_w over patches #' num2str([numbers{:}])])
        % set(gcf,'Name',['Local solution H_w over patches #' num2str([numbers{:}])])
    elseif ischarin('displ',varargin)
        i = getcharin('displ',varargin);
        figure('Name',['Local solution w_' num2str(i) ' over patches #' num2str([numbers{:}])])
        % set(gcf,'Name',['Local solution w_' num2str(i) ' over patches #' num2str([numbers{:}])])
    elseif ischarin('rotation',varargin)
        i = getcharin('rotation',varargin);
        figure('Name',['Local solution rot_w_' num2str(i) ' over patches #' num2str([numbers{:}])])
        % set(gcf,'Name',['Local solution rot_w_' num2str(i) ' over patches #' num2str([numbers{:}])])
    else
        figure('Name',['Local solution w over patches #' num2str([numbers{:}])])
        % set(gcf,'Name',['Local solution w over patches #' num2str([numbers{:}])])
    end
    clf
    n = numel(patches);
    for k=1:n
        patch = patches.patches{k};
        plot_sol(patch.S,w{patch.number},varargin{:});
        colormap(p.Results.colormap)
        if p.Results.colorbar
            colorbar
        end
        set(gca,'FontSize',p.Results.FontSize)
    end
elseif isa(patches,'Patch')
    patch = patches;
    if ischarin('sigma',varargin)
        i = getcharin('sigma',varargin);
        figure('Name',['Local solution sig_w_' num2str(i) ' over patch #' num2str(patch.number)])
        % set(gcf,'Name',['Local solution sig_w_' num2str(i) ' over patch #' num2str(patch.number)])
    elseif ischarin('epsilon',varargin)
        i = getcharin('epsilon',varargin);
        figure('Name',['Local solution eps_w_' num2str(i) ' over patch #' num2str(patch.number)])
        % set(gcf,'Name',['Local solution eps_w_' num2str(i) ' over patch #' num2str(patch.number)])
    elseif ischarin('energyint',varargin)
        figure('Name',['Local solution H_w over patch #' num2str(patch.number)])
        % set(gcf,'Name',['Local solution H_w over patch #' num2str(patch.number)])
    elseif ischarin('displ',varargin)
        i = getcharin('displ',varargin);
        figure('Name',['Local solution w_' num2str(i) ' over patch #' num2str(patch.number)])
        % set(gcf,'Name',['Local solution w_' num2str(i) ' over patch #' num2str(patch.number)])
    elseif ischarin('rotation',varargin)
        i = getcharin('rotation',varargin);
        figure('Name',['Local solution rot_w_' num2str(i) ' over patch #' num2str(patch.number)])
        % set(gcf,'Name',['Local solution rot_w_' num2str(i) ' over patch #' num2str(patch.number)])
    else
        figure('Name',['Local solution w over patch #' num2str(patch.number)])
        % set(gcf,'Name',['Local solution w over patch #' num2str(patch.number)])
    end
    clf
    plot_sol(patch.S,w{patch.number},varargin{:});
    colormap(p.Results.colormap)
    if p.Results.colorbar
        colorbar
    end
    set(gca,'FontSize',p.Results.FontSize)
end

end
