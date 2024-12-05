function myprintpdf(pathname,filename,resolution,renderer)
% function myprintpdf(pathname,filename,resolution,renderer)
% prints current figure in PDF format without tons of white space
% pathname : char designant le repertoire de sauvegarde
% filename : char designant le nom du fichier (ne pas preciser l'extension)
% resolution : double designant le nombre de dpi
% renderer : char designant le nom du rendu ('painters', 'zbuffer' or 'OpenGL')

if nargin<=2 || isempty(resolution)
    resolution = [];
end
if nargin<=3 || isempty(renderer)
    renderer = [];
end

newfile = fullfile(pathname,filename);

err = dos(['cd ' pathname]);
if err~=0
    error('the directory does not exist')
end

% The width and height of the figure are found
% The paper is set to be the same width and height as the figure
% The figure's bottom left corner is lined up with
% the paper's bottom left corner

% Set figure and paper to use the same unit
set(gcf, 'Units', 'centimeters')
set(gcf, 'PaperUnits','centimeters')

% Position of figure is of form [left bottom width height]
% We only care about width and height
pos = get(gcf,'Position');

% Set paper size to be same as figure size
set(gcf, 'PaperSize', [pos(3) pos(4)])

% Set figure to start at bottom left of paper
% This ensures that figure and paper will match up in size
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperPosition', [0 0 pos(3) pos(4)])

if ~isempty(renderer)
    set(gcf, 'Renderer',renderer);
end

if ~isempty(resolution)
    print(gcf,'-dpdf',['-r' num2str(resolution)],newfile);
else
    print(gcf,'-dpdf',newfile)
end

end
