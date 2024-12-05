function mov = evolSobolIndices(S,T,ut,alpha,varargin)
% function mov = evolSobolIndices(S,T,ut,alpha,varargin)
% Display evolution of the Sobol indices of solution ut 
% associated with the group of variables alpha in {1,..,d} with d = ndims(ut)
% S: MODEL
% T: TIMEMODEL
% ut: FunctionalBasisArray of solution u
% alpha: 1-by-s array of integers or 1-by-d logical
% - if alpha is an array of integers, indices with respect
% to variables alpha
% - if alpha is logical, indices with respect
% to variables find(alpha)
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
    figure('Name',['Sobol index of u_' num2str(i) ' for random variables #' num2str(alpha)])
    % set(gcf,'Name',['Sobol index of u_' num2str(i) ' for random variables #' num2str(alpha)])
else
    figure('Name',['Sobol index of u for random variables #' num2str(alpha)])
    % set(gcf,'Name',['Sobol index of u for random variables #' num2str(alpha)])
end
clf
set(gcf,'color','w')

sz = ut.sz;
d = ndims(ut);
sut = SensitivityAnalysis.sobolIndices(ut,alpha,d);
sut = reshape(sut,sz);
sut = TIMEMATRIX(sut,T);
sut = unfreevector(S,sut)-calc_init_dirichlet(S)*one(T);

sut = setevolparam(sut,'colormap',p.Results.colormap,'colorbar',p.Results.colorbar,...
    'view',p.Results.view,'camup',p.Results.camup,'FontSize',p.Results.FontSize);
frame = evol_sol(sut,S,'rescale',p.Results.rescale,varargin{:}); % save the frames

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
