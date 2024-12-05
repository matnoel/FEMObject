function mysaveas(pathname,filename,formats,renderer)
% function mysaveas(pathname,filename,formats,renderer)
% save current figure
% pathname : char containing the path name
% filename : char containing the file name (without extension)
% formats : cell array, each cell containing the specified format
% renderer : char containing the renderer ('painters', 'zbuffer' or 'opengl')

if nargin==2 || isempty(formats)
    formats = {'epsc'};
elseif nargin>=3 && isa(formats,'char')
    formats = {formats};
end
if nargin<=3 || isempty(renderer)
    renderer = [];
end

newfile = fullfile(pathname,filename);

err = dos(['cd ' pathname]);
if err~=0
    error('the directory does not exist')
end

set(gcf,'PaperPositionMode','auto')
if ~isempty(renderer)
    set(gcf,'Renderer',renderer)
end

for k=1:length(formats)
    saveas(gcf,newfile,[formats{k}])
end

if ischarin('jpeg',formats) && ~ischarin('eps',formats) && ~ischarin('eps2',formats) && ~ischarin('epsc',formats) && ~ischarin('epsc2',formats)
    com = ['convert ' newfile '.jpg -compress jpeg eps2:' newfile '.eps &'];
    dos(com);
end

if ischarin('pdf',formats)
    com = ['pdfcrop ' newfile '.pdf &'];
    dos(com);
end

end
