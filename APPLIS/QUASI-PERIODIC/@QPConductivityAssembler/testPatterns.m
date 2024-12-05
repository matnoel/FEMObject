function Assembler = testPatterns(varargin)

% Hard-coded bounds
maxValue = 100 ;
minValue = 1 ;
valGen = @() minValue-1+randi(maxValue-minValue+1) ;

% Create patterns (random parameters)
% uniform
patterns(1) = struct('name','uniform','value',valGen(),...
    'size',[],'center',[],'offset',[]);
% disc
patterns(2) = struct('name','disc','value',valGen(),...
    'size',rand(1),'center',rand(1,2),'offset',[]) ;
% rectangle'value',valGen(),
patterns(3) = struct('name','rectangle','value',valGen(),...
    'size',rand(1,2),'center',rand(1,2),'offset',[]) ;
% chevron
patterns(4) = struct('name','chevron','value',valGen(),...
    'size',circshift([1 rand(1)],[0 randi(2)]),'center',rand(1,2),'offset',rand(1)) ;
% bar
patterns(5) = struct('name','bar','value',valGen(),...
    'size',circshift([1 rand(1)],[0 randi(2)]),'center',rand(1,2),'offset',[]) ;
% cross
patterns(6) = struct('name','cross','value',valGen(),...
    'size',rand(1,2),'center',rand(1,2),'offset',rand(1,2)) ;

% Pattern table and distribution
patternNb = numel(patterns) ;
patternTable = eye(patternNb) ; % one pattern -> one field
patternTable(1,:) = 1 ; % add 'uniform' pattern to every field
distributor = @(n) num2cell((1:n)') ;
cellNum = subplot_format(patternNb) ;

% Create assembler
Assembler = QPConductivityAssembler.createRandom('patterns',patterns,...
    'patternTable',patternTable,'distributor',distributor,...
    'cellNum',cellNum,varargin{:}) ;
end