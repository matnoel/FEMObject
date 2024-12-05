function pb = updateUpdater(pb,lhs,rhs)
% pb = updateUpdater(pb,lhs,rhs)

%TODO: update stagnation from tolerance ?

if nargin < 3
    rhs = getRHSOperator(pb) ;
    if nargin < 2 || isempty(lhs)
        lhs = getLHSOperator(pb) ;
    end
end

updater = getUpdater(pb) ;
if isempty(updater)
    updater = defaultUpdater(pb) ;
end

[lhs,rhs] = convertTensors(lhs,rhs) ;
updater.A = lhs ;
updater.b = rhs ;
updater.tolerance = getTolerance(pb) ;
stagTol = min(1e-6,updater.tolerance/100) ;
updater.stagnation = stagTol ;
updater.stagnationCore = stagTol ;
updater.stagnationTSpace = stagTol ;
updater.display = getVerbose(pb) ;

pb = setUpdater(pb,updater) ;
end