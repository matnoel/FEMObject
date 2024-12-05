function K = getCompleteK(pb,fullModel)
% K = getCompleteK(pb,fullModel)
% Wrapper used by solveFE in case of QPDDProblem.

if nargin < 2
    fullModel = [] ;
end

K = assembleGlobalLocal(pb,pb.K,{pb.patch(:).K},fullModel) ;
    
end