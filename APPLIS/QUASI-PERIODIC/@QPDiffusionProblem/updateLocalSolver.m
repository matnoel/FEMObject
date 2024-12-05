function pb = updateLocalSolver(pb,lhs,rhs)
% pb = updateLocalSolver(pb,lhs,rhs)

%TODO: update stagnation from tolerance ?

if nargin < 3
    rhs = getRHSOperator(pb) ;
    if nargin < 2 || isempty(lhs)
        lhs = getLHSOperator(pb) ;
    end
end

localSolver = getLocalSolver(pb) ;
if isempty(localSolver)
    localSolver = defaultLocalSolver(pb) ;
end

[lhs,rhs] = convertTensors(lhs,rhs) ;
localSolver.A = lhs ;
localSolver.b = rhs ;

pb = setLocalSolver(pb,localSolver) ;
end