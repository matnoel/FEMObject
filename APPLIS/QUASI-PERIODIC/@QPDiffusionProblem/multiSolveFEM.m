function [solutions,outputs] = multiSolveFEM(DP,distributions,mesherMethod)
% [solutions,outputs] = multiSolveFEM(DP,distributions,mesherMethod)

nbDist = size(distributions,1) ;
solutions = cell(nbDist,1) ;

for iter = 1:nbDist
    [currentPb,assTime] = updateDistribution(DP,distributions(iter,:)) ;
    [solutions{iter},outputs(iter)] = solveFEM(currentPb,mesherMethod) ;
    % Add conductivity assembling time to operator's
    outputs(iter).time(2) = outputs(iter).time(2) + assTime ;
end

end