function [solution,output,pb] = solve(pb,orders2assemble)
% [solution,output,pb] = solve(pb,orders2assemble)

if nargin == 1
    orders2assemble = [] ;
end

initTime = 0 ;

% Check whether assembling is required
patches = getPatches(pb) ;
requirePatcheAssembling = ~isempty(patches) && ...
    (isempty(getLHSOperator(patches)) || isempty(getRHSOperator(patches)) ) ;
requireAssembling = ~isempty(orders2assemble) || requirePatcheAssembling...
    || isempty(getLHSOperator(pb)) || isempty(getRHSOperator(pb)) ;

if requireAssembling
    if isempty(orders2assemble)
        warning('Forcing assembling because required operators are missing.')
        orders2assemble = 1:getOrder(pb) ;
    end
    [pb,assTime] = assembleOperators(pb,orders2assemble) ;
    initTime = initTime + assTime ;
end

% Solve according to whether there are patches
if isempty(getPatches(pb))
    % Compression
    tic ;
    compressor = Truncator('tolerance',getTolSVD(pb)) ;
    lhs = truncate(compressor,getLHSOperator(pb)) ;
    rhs = truncate(compressor,getRHSOperator(pb)) ;
    pb = updateGreedySolver(pb,lhs,rhs) ;
    initTime = initTime + toc ;
    % Resolution
    greedy = getGreedySolver(pb) ;
    [solution,output] = solve(greedy) ;
    output.initTime = initTime ;
else
    [solution,output,pb] = solveDD(pb) ;
    output.initTime = output.initTime + initTime ;
end

end