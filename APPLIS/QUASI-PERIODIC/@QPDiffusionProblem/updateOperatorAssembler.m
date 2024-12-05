function pb = updateOperatorAssembler(pb,operatorAssembler)
% pb = updateOperatorAssembler(pb,operatorAssembler)

if nargin == 1
    operatorAssembler = getOperatorAssembler(pb) ;
end

pb = setOperatorAssembler(pb,operatorAssembler) ;
pb = updateGreedySolver(pb) ;

end