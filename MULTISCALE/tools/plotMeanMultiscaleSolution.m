function plotMeanMultiscaleSolution(glob,patches,interfaces,U,w,varargin)
% function plotMeanMultiscaleSolution(glob,patches,interfaces,U,w,varargin)
% Display the mean (mathematical expectation) of multiscale solution u=(U,w)
% glob: Global
% patches: Patches
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
    figure('Name',['Mean of u_' num2str(i) '=(U_' num2str(i) ',w_' num2str(i) ') over domain'])
    % set(gcf,'Name',['Mean of u_' num2str(i) '=(U_' num2str(i) ',w_' num2str(i) ') over domain'])
else
    figure('Name','Mean of u=(U,w) over domain')
    % set(gcf,'Name','Mean of u=(U,w) over domain')
end
clf

if U.sz(1)==glob.S.nbddl
    P_out = calcProjection(glob,'free',false);
else
    P_out = glob.P_out;
end

plot_sol(glob.S_out,mean(U*P_out')',varargin{:});
for k=1:n
    patch = patches.patches{k};
    interface = interfaces.interfaces{k};
    plot_sol(patch.S,mean(w{patch.number})',varargin{:});
    if ~ischarin('sigma',varargin) && ~ischarin('epsilon',varargin) && ~ischarin('energyint',varargin)
        plot_sol(interface.S,mean(w{patch.number}*interface.P_patch')','FaceColor','none','EdgeColor','k',varargin{:});
    end
end
colormap(p.Results.colormap)
if p.Results.colorbar
    colorbar
end
set(gca,'FontSize',p.Results.FontSize)

end
