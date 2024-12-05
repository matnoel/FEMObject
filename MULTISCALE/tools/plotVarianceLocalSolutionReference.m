function plotVarianceLocalSolutionReference(patches,w_ref,varargin)
% function plotVarianceLocalSolutionReference(patches,w_ref,varargin)
% Display the variance of reference local solution w_ref
% patches: Patches or Patch
% w_ref: FunctionalBasisArray of reference local solution w

p = ImprovedInputParser;
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

varargin = delcharin({'colorbar','colormap','FontSize'},varargin);

if isa(patches,'Patches')
    numbers = getnumber(patches);
    if ischarin('displ',varargin)
        i = getcharin('displ',varargin);
        figure('Name',['Variance of w_ref_' num2str(i) ' over patches #' num2str([numbers{:}])])
        % set(gcf,'Name',['Variance of w_ref_' num2str(i) ' over patches #' num2str([numbers{:}])])
    else
        figure('Name',['Variance of w_ref over patches #' num2str([numbers{:}])])
        % set(gcf,'Name',['Variance of w_ref over patches #' num2str([numbers{:}])])
    end
    clf
    n = numel(patches);
    for k=1:n
        patch = patches.patches{k};
        plot_sol(patch.S,variance(w_ref{patch.number})',varargin{:});
        colormap(p.Results.colormap)
        if p.Results.colorbar
            colorbar
        end
        set(gca,'FontSize',p.Results.FontSize)
    end
elseif isa(patches,'Patch')
    patch = patches;
    if ischarin('displ',varargin)
        i = getcharin('displ',varargin);
        figure('Name',['Variance of w_ref_' num2str(i) ' over patch #' num2str(patch.number)])
        % set(gcf,'Name',['Variance of w_ref_' num2str(i) ' over patch #' num2str(patch.number)])
    else
        figure('Name',['Variance of w_ref over patch #' num2str(patch.number)])
        % set(gcf,'Name',['Variance of w_ref over patch #' num2str(patch.number)])
    end
    clf
    plot_sol(patch.S,variance(w_ref{patch.number})',varargin{:});
    colormap(p.Results.colormap)
    if p.Results.colorbar
        colorbar
    end
    set(gca,'FontSize',p.Results.FontSize)
end

end
