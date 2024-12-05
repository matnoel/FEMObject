function plotSensitivityIndicesLagrangeMultiplier(interfaces,lambda,alpha,varargin)
% function plotSensitivityIndicesLagrangeMultiplier(interfaces,lambda,alpha,varargin)
% Display the sensitivity indices of Lagrange multiplier lambda associated with the
% group of variables alpha in {1,..,d} with d = ndims(lambda)
% with respect to the single random variable in dimension dim
% interfaces: Interfaces or Interface
% lambda: FunctionalBasisArray of Lagrange multiplier lambda
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

if isa(interfaces,'Interfaces')
    numbers = getnumber(interfaces);
    if ischarin('displ',varargin)
        i = getcharin('displ',varargin);
        figure('Name',['Sensitivity index of lambda_' num2str(i) ' over interfaces #' num2str([numbers{:}]) ' for random variables #' num2str(alpha)])
        % set(gcf,'Name',['Sensitivity index of lambda_' num2str(i) ' over interfaces #' num2str([numbers{:}]) ' for random variables #' num2str(alpha)])
    else
        figure('Name',['Sensitivity index of lambda over interfaces #' num2str([numbers{:}]) ' for random variables #' num2str(alpha)])
        % set(gcf,'Name',['Sensitivity index of lambda over interfaces #' num2str([numbers{:}]) ' for random variables #' num2str(alpha)])
    end
    clf
    n = numel(interfaces);
    for k=1:n
        interface = interfaces.interfaces{k};
        v = variance(lambda{interface.number});
        s = varianceConditionalExpectation(lambda{interface.number},alpha)'./max(v);
        plot_sol(interface.S,s,'EdgeColor','interp',varargin{:});
        colormap(p.Results.colormap)
        if p.Results.colorbar
            colorbar
        end
        set(gca,'FontSize',p.Results.FontSize)
    end
elseif isa(interfaces,'Interface')
    interface = interfaces;
    if ischarin('displ',varargin)
        i = getcharin('displ',varargin);
        figure('Name',['Sensitivity index of lambda_' num2str(i) ' over interface #' num2str(interface.number) ' for random variables #' num2str(alpha)])
        % set(gcf,'Name',['Sensitivity index of lambda_' num2str(i) ' over interface #' num2str(interface.number) ' for random variables #' num2str(alpha)])
    else
        figure('Name',['Sensitivity index of lambda over interface #' num2str(interface.number) ' for random variables #' num2str(alpha)])
        % set(gcf,'Name',['Sensitivity index of lambda over interface #' num2str(interface.number) ' for random variables #' num2str(alpha)])
    end
    clf
    v = variance(lambda{interface.number});
    s = varianceConditionalExpectation(lambda{interface.number},alpha)'./max(v);
    plot_sol(interface.S,s,'EdgeColor','interp',varargin{:});
    colormap(p.Results.colormap)
    if p.Results.colorbar
        colorbar
    end
    set(gca,'FontSize',p.Results.FontSize)
end

end
