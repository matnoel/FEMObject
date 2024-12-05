function pb = updateConductivity(pb,K,distribution)
% pb = updateConductivity(pb,K,distribution)
% Distribution is optional and helps to compute conductivity bounds.

if nargin < 3
    distribution = [] ;
if nargin == 1
    K = getConductivity(pb) ;
end
end

cA = getConductivityAssembler(getOperatorAssembler(pb)) ;
cA = setConductivity(cA,K) ;
if ~isempty(distribution)
    cA = setDistribution(cA,distribution) ;
end
pb = updateConductivityAssembler(pb,cA) ;

end