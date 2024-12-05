function pb = updateGreedySolver(pb,lhs,rhs)
% pb = updateGreedySolver(pb,lhs,rhs)

if nargin < 3
    rhs = getRHSOperator(pb) ;
    if nargin < 2 || isempty(lhs)
        lhs = getLHSOperator(pb) ;
    end
end

% Update greedy

greedy = getGreedySolver(pb) ;
if isempty(greedy)
    greedy = defaultGreedySolver(pb) ;
end

[lhs,rhs] = convertTensors(lhs,rhs) ;
greedy.A = lhs ;
greedy.b = rhs ;
greedy.tolerance = getTolerance(pb) ;

% Updater
pb = updateUpdater(pb,lhs,rhs) ;
up = getUpdater(pb) ;
greedy.update = @(x) updateTSpaceByALS(up,x) ;

pb = setGreedySolver(pb,greedy) ;

% Update dependencies
pb = updateLocalSolver(pb,lhs,rhs) ;
end