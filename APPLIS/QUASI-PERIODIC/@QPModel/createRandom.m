function model = createRandom(varargin)
% model = createRandom(varargin)

% Hard coded bounds
rgCSz = [1e-1 5] ; % range for cell size
rgCNu = [1 4] ; % range for cell numbers along axes
rgElNb = [3 30] ; % range for elements number

% Parse varargin and and random values to missing properties
if ~ischarin('order',varargin)
    varargin = setcharin('order',varargin,1+randi(2)) ;
end
if ~ischarin('cellNum',varargin)
    varargin = setcharin('cellNum',varargin,rgCNu(1)-1+randi(rgCNu(2),1,2)) ;
end
if ~ischarin('cellSize',varargin)
    cellSize = max(rgCSz(2)*rand(1,2),rgCSz(1)) ;
    varargin = setcharin('cellSize',varargin,cellSize) ;
end
if ~ischarin('elementSize',varargin)
    elementNumber = max( randi(rgElNb(2),1,2) ,rgElNb(1)) ;
    elementSize = cellSize./elementNumber ;
    varargin = setcharin('elementSize',varargin,elementSize) ;
end
if ~ischarin('verbose',varargin)
    varargin = setcharin('verbose',varargin,false) ;
end

model = QPModel(varargin{:}) ;
end