function DP = createRandom(varargin)
% DP = createRandom(varargin)

% Parse varargin and and random values to missing properties
if ~ischarin('operatorAssembler',varargin)
    rdOpAss = QPDiffusionAssembler.createRandom(varargin{:}) ;
    varargin = setcharin('operatorAssembler',varargin,rdOpAss) ;
end
if ~ischarin('verbose',varargin)
    varargin = setcharin('verbose',varargin,false) ;
end

DP = QPDiffusionProblem(varargin{:}) ;
end