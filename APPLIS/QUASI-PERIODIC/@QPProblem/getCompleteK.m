function K = getCompleteK(pb,fullModel)
% K = getCompleteK(pb,fullModel)
% Wrapper used by solveFE in case of QPDDProblem.

if nargin < 2
    fullModel = [] ;
end

K = pb.K ;
if isa(K,'TuckerLikeTensor')
    
    K = untensorize(pb.model,K) ;
    qpCoord = smooth(pb.model,getDomainCoord(pb.model)) ;
    K = relativeSort(K,qpCoord,getcoord(getnode(fullModel)),...
        getfemobjectoptions('tolerancepoint')) ;
    
elseif isnumeric(K)
    
    if isscalar(K)
        K = K(ones(getnbnode(fullModel),1)) ;
    else
        K = smooth(pb.model,K) ;
        qpCoord = getDiscon2Con(pb.model)*getDomainCoord(pb.model) ;
        K = relativeSort(K,qpCoord,getcoord(getnode(fullModel)),...
            getfemobjectoptions('tolerancepoint')) ;
    end
    
end
end