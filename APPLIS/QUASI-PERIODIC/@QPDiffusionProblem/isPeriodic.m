function bool = isPeriodic(pb,direction)
% bool = isPeriodic(pb,direction)
% Calls to QPDiffusionAssembler homonym method.

if nargin < 2
    direction = [1 2] ;
end

bool = isPeriodic(getOperatorAssembler(pb),1:2) ;

end