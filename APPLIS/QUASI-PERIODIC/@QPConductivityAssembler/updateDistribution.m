function [Assembler,assemblerTime] = updateDistribution(Assembler,newDistribution)
% [Assembler,assemblerTime] = updateDistribution(Assembler,newDistribution)

if nargin == 1 || isempty(newDistribution)
    newDistribution = distribute(Assembler) ;
end

Assembler = setDistribution(Assembler,newDistribution) ;
[Assembler,assemblerTime] = assemble(Assembler) ;

end