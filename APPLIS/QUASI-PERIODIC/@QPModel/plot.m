function [] = plot(model,t,varargin)
% [] = plot(model,t,varargin)

smoothRequested = ~ischarin('nosmooth',varargin) && ~ischarin('rough',varargin) ;

coord = getcharin('coord',varargin,getDomainCoord(model)) ;
% Process according to input's class
if isa(t,'TuckerLikeTensor')
    t = untensorize(model,t,smoothRequested) ;
elseif iscell(t) % cell recursion
    subplotFormat = subplot_format(numel(t)) ;
    for i = 1:numel(t)
        subplot(subplotFormat(1),subplotFormat(2),i) ;
        plot(model,t{i},varargin{:}) ;
    end
    return
end

if smoothRequested
    coord = smooth(model,coord) ; % smooth has built-in safety to
    t = smooth(model,t) ;         % recognize already smooth fields
end
xyz = {full(coord(:,1)) ; full(coord(:,2)) ; full(t)} ;

% Process varargin options
varargin = delonlycharin('nosmooth',varargin) ; % remove before passing it
varargin = delonlycharin('rough',varargin) ;    % to trisurf
varargin = delcharin('coord',varargin) ;
if ~ischarin('EdgeColor',varargin)
    varargin = setcharin('EdgeColor',varargin,'interp') ;
end

if ~ischarin('FaceColor',varargin)
    varargin = setcharin('FaceColor',varargin,'interp') ;
end

if ischarin('cellBorder',varargin)
    drawCellBorder = true ;
    varargin = delonlycharin('cellBorder',varargin) ;
else
    drawCellBorder = false ;
end

% Draw plot
gridDelaunay = delaunay(xyz{1},xyz{2}) ;
trisurf(gridDelaunay,xyz{1},xyz{2},xyz{3},varargin{:}) ;
axis([min(xyz{1}) max(xyz{1}) min(xyz{2}) max(xyz{2})]) ;
% Previous line restrict axis and set view from above, as view([0 0 1])
% To keep default view, use
% axis([min(xyz{1}) max(xyz{1}) min(xyz{2}) max(xyz{2}) min(xyz{3}) max(xyz{4})]) ;

% Draw cells borders
if drawCellBorder
    cellNum = getCellNum(model) ;
    domainSize = max(coord);
    cellSize = domainSize./cellNum;
    
    border = cell(2,1) ;
    for d = 1:2
        nd = 3-d ; % not direction (the other one)
        borderCoord = (1:cellNum(d)-1)'*cellSize(d) ;
        borderCoord = kron(borderCoord,[1;1]) ; % replicte each row twice
        borderWidth = [0 ; domainSize(nd)] ;
        borderWidth = kron(ones(cellNum(d)-1,1),borderWidth) ; % replicate
        % whole column cellNum(d) times
        border{d} = [borderCoord borderWidth] ;
    end
    border = [border{1} ; circshift(border{2},[0 1])] ; % concat with permutation
    borderHeight = (max(xyz{3})+1)*ones(size(border,1),1) ;
    border = [border borderHeight] ;
    
    hold on
    for i=1:size(border,1)/2
        j = [2*i-1 ; 2*i] ;
        plot3(border(j,1),border(j,2),border(j,3),'--w','LineWidth',2)
    end
    hold off
end

end