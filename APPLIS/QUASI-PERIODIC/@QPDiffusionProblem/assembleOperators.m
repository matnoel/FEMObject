function [DP,assemblerTime] = assembleOperators(DP,orders2assemble)
% [DP,assemblerTime] = assembleOperators(DP,orders2assemble)

if nargin == 1
    orders2assemble = 1:getOrder(DP) ;
end

opAss = getOperatorAssembler(DP) ;
[opAss,assemblerTime] = assemble(opAss,orders2assemble) ;
DP = updateOperatorAssembler(DP,opAss) ;
end