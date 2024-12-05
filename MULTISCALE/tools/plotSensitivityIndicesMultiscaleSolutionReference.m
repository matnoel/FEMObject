function plotSensitivityIndicesMultiscaleSolutionReference(glob,patches,interfaces,U_ref,w_ref,alpha,varargin)
% function plotSensitivityIndicesMultiscaleSolutionReference(glob,patches,interfaces,U_ref,w_ref,alpha,varargin)
% Display the sensitivity indices of reference multiscale solution u_ref=(U_ref,w_ref) 
% associated with the group of variables alpha in {1,..,d} with d = ndims(u_ref)
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
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

varargin = delcharin({'colorbar','colormap','FontSize'},varargin);
n = numel(patches);

if ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Sensitivity index of u_ref_' num2str(i) '=(U_ref_' num2str(i) ',w_ref_' num2str(i) ') over domain for random variables #' num2str(alpha)])
    % set(gcf,'Name',['Sensitivity index of u_ref_' num2str(i) '=(U_ref_' num2str(i) ',w_ref_' num2str(i) ') over domain for random variables #' num2str(alpha)])
else
    figure('Name',['Sensitivity index of u_ref=(U_ref,w_ref) over domain for random variables #' num2str(alpha)])
    % set(gcf,'Name',['Sensitivity index of u_ref=(U_ref,w_ref) over domain for random variables #' num2str(alpha)])
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

sU_ref = varianceConditionalExpectation(U_ref,alpha)'./max_var;
sU_ref = unfreevector(S_out,sU_ref)-calc_init_dirichlet(S_out);
plot_sol(S_out,sU_ref,varargin{:});
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
    colorbar
end
set(gca,'FontSize',p.Results.FontSize)

end
