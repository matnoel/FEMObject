function pb = updateModel(pb,model)
% DP = updateModel(DP,model)
if nargin == 1
    model = getModel(pb) ;
end

operatorAssembler = updateModel(getOperatorAssembler(pb),model) ;
pb = updateOperatorAssembler(pb,operatorAssembler) ;

% update maxIterations
greedy = getGreedySolver(pb) ;
greedy.maxIterations = max(getCellNb(model),greedy.maxIterations) ;
pb = setGreedySolver(pb,greedy) ;
end