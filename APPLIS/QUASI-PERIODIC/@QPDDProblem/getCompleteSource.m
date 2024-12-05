function source = getCompleteSource(pb,fullModel)
% source = getCompleteSource(pb,fullModel)
% Wrapper used by solveFE in case of QPDDProblem.

if nargin < 2
    fullModel = [] ;
end

source = pb.source ;
if ischar(source) % then is corrector source
    return
end

source = assembleGlobalLocal(pb,pb.source,{pb.patch(:).source},fullModel) ;
    
end