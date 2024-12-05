function [solutions,outputs,pb] = multiSolve(pb,dist)
% [solutions,outputs] = multiSolve(DP,dist)

if iscell(dist)
    nb = numel(dist) ;
    randomize = false ;
elseif isnumeric(dist)
    randomize = true ;
    nb = dist ;
else
    error('Wrong type for second input argument.')
end

pbOrder = getOrder(pb) ;
solutions = cell(nb,1) ;
outputs = QPDiffusionProblem.qpOutputStructure([nb 1]) ;

% Initial run
if randomize
    pb = updateDistribution(pb) ;
else
    pb = updateDistribution(pb,dist{1}) ;
end
[solutions{1},outputs(1),pb] = solve(pb,1:pbOrder) ;

% Iterative runs
for iter = 2:nb
    pb = setInitialPoint(pb,solutions{iter-1}) ;
    if randomize
        [pb,assTime] = updateDistribution(pb) ;
    else
        [pb,assTime] = updateDistribution(pb,dist{iter}) ;
    end
    [solutions{iter},outputs(iter),pb] = solve(pb,1:pbOrder-1) ;
    % Add conductivity assembling time to operator's
    outputs(iter).initTime = outputs(iter).initTime + assTime ;
end

end