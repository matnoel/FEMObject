function [DP,assemblerTime] = updateDistribution(DP,distribution)
% [DP,assemblerTime] = updateDistribution(DP,distribution)

if nargin==1
    distribution = [] ; % QPConductivityAssembler will know what to do
end

KAss = getConductivityAssembler(DP) ;
[KAss,assemblerTime] = updateDistribution(KAss,distribution) ;
DP = updateConductivityAssembler(DP,KAss) ;

end