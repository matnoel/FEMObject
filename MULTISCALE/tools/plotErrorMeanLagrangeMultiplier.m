function plotErrorMeanLagrangeMultiplier(interfaces,lambda,lambda_ref,varargin)
% function plotErrorMeanLagrangeMultiplier(interfaces,lambda,lambda_ref,varargin)
% Display the relative error in the mean (mathematical expectation) of Lagrange multiplier lambda
% with respect to reference Lagrange multiplier lambda_ref
% interfaces: Interfaces
% lambda: FunctionalBasisArray of Lagrange multiplier lambda
% lambda_ref: FunctionalBasisArray of reference Lagrange multiplier lambda

p = ImprovedInputParser;
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

varargin = delcharin({'colorbar','colormap','FontSize'},varargin);
n = numel(interfaces);

if ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Relative error in mean of lambda_' num2str(i) ' w.r.t. lambda_ref_' num2str(i)])
    % set(gcf,'Name',['Relative error in mean of lambda_' num2str(i) ' w.r.t. lambda_ref_' num2str(i)])
else
    figure('Name','Relative error in mean of lambda w.r.t. lambda_ref')
    % set(gcf,'Name','Relative error in mean of lambda w.r.t. lambda_ref')
end
clf

for k=1:n
    interface = interfaces.interfaces{k};
    err_interface = abs(mean(lambda{interface.number})-mean(lambda_ref{interface.number}))'/max(mean(lambda{interface.number}));
    plot_sol(interface.S,err_interface,'EdgeColor','interp',varargin{:});
    colormap(p.Results.colormap)
    if p.Results.colorbar
        colorbar
    end
    set(gca,'FontSize',p.Results.FontSize)
end

end
