function mov = evolMeanGlobalSolution(glob,T,Ut,varargin)
% function mov = evolMeanGlobalSolution(glob,T,Ut,varargin)
% Display evolution of the mean (mathematical expectation) of global solution Ut
% glob: Global
% T: TIMEMODEL
% Ut: FunctionalBasisArray of global solution U
% mov: movie

p = ImprovedInputParser;
addParameter(p,'rescale',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colorbar',true,@(x) islogical(x) || ischar(x));
addParameter(p,'colormap','default',@(x) isnumeric(x) || ischar(x));
addParameter(p,'view',[],@isnumeric);
addParameter(p,'camup','auto',@(x) isnumeric(x) || ischar(x));
addParameter(p,'FontSize',16,@isscalar);
addParameter(p,'filename','solution',@ischar);
addParameter(p,'pathname','./',@ischar);
addParameter(p,'formats',{'avi','mp4'},@(x) ischar(x) || iscell(x));
addParameter(p,'FrameRate',30,@isnumeric);
addParameter(p,'Quality',100,@isnumeric);
parse(p,varargin{:})

varargin = delcharin({'rescale','colorbar','colormap','view','camup','FontSize',...
    'filename','pathname','formats','FrameRate','Quality'},varargin);
if isa(p.Results.formats,'char')
    p.Results.formats = {p.Results.formats};
end

% Create a new figure
if ischarin('displ',varargin)
    i = getcharin('displ',varargin);
    figure('Name',['Mean of U_' num2str(i) ' over fictitious domain'])
    % set(gcf,'Name',['Mean of U_' num2str(i) ' over fictitious domain'])
else
    figure('Name','Mean of U over fictitious domain')
    % set(gcf,'Name','Mean of U over fictitious domain')
end
clf
set(gcf,'color','w')

sz = Ut.sz;
Ut = mean(Ut);
Ut = reshape(Ut,sz);
Ut = TIMEMATRIX(Ut,T);

Ut = setevolparam(Ut,'colormap',p.Results.colormap,'colorbar',p.Results.colorbar,...
    'view',p.Results.view,'camup',p.Results.camup,'FontSize',p.Results.FontSize);
frame = evol_sol(Ut,glob.S,'rescale',p.Results.rescale,varargin{:}); % save the frames

% Create movie file
mov = cell(1,length(p.Results.formats));
for i=1:length(p.Results.formats)
    if strcmp(p.Results.formats{i},'avi')
        mov{i} = VideoWriter(fullfile(p.Results.pathname,p.Results.filename));
    elseif strcmp(p.Results.formats{i},'mp4')
        mov{i} = VideoWriter(fullfile(p.Results.pathname,p.Results.filename),'MPEG-4');
    elseif strcmp(p.Results.formats{i},'mj2')
        mov{i} = VideoWriter(fullfile(p.Results.pathname,p.Results.filename),'Motion JPEG 2000');
    end
    mov{i}.FrameRate = p.Results.FrameRate;
    mov{i}.Quality = p.Results.Quality;
    open(mov{i});
    writeVideo(mov{i},frame); % add the frames to the movie
    close(mov{i});
end

end
