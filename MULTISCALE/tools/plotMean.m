function plotMean(S,u,varargin)
% function plotMean(S,u,varargin)
% Display the mean (mathematical expectation) of solution u
% S: MODEL
% u: FunctionalBasisArray of solution u

p = ImprovedInputParser;
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
parse(p,varargin{:})

varargin = delcharin({'colorbar','colormap','FontSize'},varargin);

if ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Mean of u_' num2str(i)])
    % set(gcf,'Name',['Mean of u_' num2str(i)])
else
    figure('Name','Mean of u')
    % set(gcf,'Name','Mean of u')
end
clf

plot_sol(S,mean(u)',varargin{:});
colormap(p.Results.colormap)
if p.Results.colorbar
    colorbar
end
set(gca,'FontSize',p.Results.FontSize)

end
