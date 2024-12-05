function plotStdLagrangeMultiplier(interfaces,lambda,varargin)
% function plotStdLagrangeMultiplier(interfaces,lambda,varargin)
% Display the standard deviation of Lagrange multiplier lambda
% interfaces: Interfaces or Interface
% lambda: FunctionalBasisArray of Lagrange multiplier lambda

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
        figure('Name',['Standard deviation of lambda_' num2str(i) ' over interfaces #' num2str([numbers{:}])])
        % set(gcf,'Name',['Standard deviation of lambda_' num2str(i) ' over interfaces #' num2str([numbers{:}])])
    else
        figure('Name',['Standard deviation of lambda over interfaces #' num2str([numbers{:}])])
        % set(gcf,'Name',['Standard deviation of lambda over interfaces #' num2str([numbers{:}])])
    end
    clf
    n = numel(interfaces);
    for k=1:n
        interface = interfaces.interfaces{k};
        plot_sol(interface.S,std(lambda{interface.number})','EdgeColor','interp',varargin{:});
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
        figure('Name',['Standard deviation of lambda_' num2str(i) ' over interface #' num2str(interface.number)])
        % set(gcf,'Name',['Standard deviation of lambda_' num2str(i) ' over interface #' num2str(interface.number)])
    else
        figure('Name',['Standard deviation of lambda over interface #' num2str(interface.number)])
        % set(gcf,'Name',['Standard deviation of lambda over interface #' num2str(interface.number)])
    end
    clf
    plot_sol(interface.S,std(lambda{interface.number})','EdgeColor','interp',varargin{:});
    colormap(p.Results.colormap)
    if p.Results.colorbar
        colorbar
    end
    set(gca,'FontSize',p.Results.FontSize)
end

end
