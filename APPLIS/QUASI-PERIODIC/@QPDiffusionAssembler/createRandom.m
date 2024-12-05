function Assembler = createRandom(varargin)
% Assembler = createRandom(varargin)

% Hard-coded bounds
maxSrcVal = 100 ;

% Parse varargin and and random values to missing properties
if ~ischarin('conductivityAssembler',varargin)
    kAss = QPConductivityAssembler.createRandom(varargin{:}) ;
    kAss = assemble(kAss) ;
    varargin = setcharin('conductivityAssembler',varargin,kAss) ;
end
if ~ischarin('source',varargin)
    model = getModel(getcharin('conductivityAssembler',varargin)) ;
    source = randomizeSource(model,maxSrcVal) ;
    varargin = setcharin('source',varargin,source) ;
end

Assembler = QPDiffusionAssembler(varargin{:}) ;
end