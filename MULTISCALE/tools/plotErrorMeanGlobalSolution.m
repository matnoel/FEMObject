function plotErrorMeanGlobalSolution(glob,U,U_ref,varargin)
% function plotErrorMeanGlobalSolution(glob,U,U_ref,varargin)
% Display the relative error in the mean (mathematical expectation) of global solution U
% with respect to reference global solution U_ref
% glob: Global
% U: FunctionalBasisArray of global solution U
% U_ref: FunctionalBasisArray of reference global solution U

p = ImprovedInputParser;
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

varargin = delcharin({'colorbar','colormap','FontSize'},varargin);

if ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Relative error in mean of U_' num2str(i) ' w.r.t. U_ref_' num2str(i)])
    % set(gcf,'Name',['Relative error in mean of U_' num2str(i) ' w.r.t. U_ref_' num2str(i)])
else
    figure('Name','Relative error in mean of U w.r.t. U_ref')
    % set(gcf,'Name','Relative error in mean of U w.r.t. U_ref')
end
clf

if U.sz(1)==glob.S.nbddl
    P_out = calcProjection(glob,'free',false);
else
    P_out = glob.P_out;
end

err_out = abs(mean(U*P_out')-mean(U_ref))'./max(mean(U_ref));
err_out = unfreevector(glob.S_out,err_out)-calc_init_dirichlet(glob.S_out);
plot_sol(glob.S_out,err_out,varargin{:});
colormap(p.Results.colormap)
if p.Results.colorbar
    colorbar
end
set(gca,'FontSize',p.Results.FontSize)

end
