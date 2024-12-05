function DP = updateConductivityAssembler(DP,conductivityAssembler)
% DP = updateConductivityAssembler(DP,conductivityAssembler)
if nargin == 1
    conductivityAssembler = getConductivityAssembler(DP) ;
end

operatorAssembler = updateConductivityAssembler(getOperatorAssembler(DP),...
    conductivityAssembler) ;
DP = updateOperatorAssembler(DP,operatorAssembler) ;
end