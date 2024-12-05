function source = getCompleteSource(pb,fullModel)
% source = getCompleteSource(pb,fullModel)
% Wrapper used by solveFE in case of QPDDProblem.

if nargin < 2
    fullModel = [] ;
end

source = pb.source ;

if isa(source,'TuckerLikeTensor')
    
    source = untensorize(pb.model,source) ;
    qpCoord = getDiscon2Con(pb.model)*getDomainCoord(pb.model) ;
    source = relativeSort(source,qpCoord,getcoord(getnode(fullModel)),...
        getfemobjectoptions('tolerancepoint')) ;
    
elseif ischar(source)
    
    return ;
    
elseif isnumeric(source)
    
    if isscalar(source)
        source = source(ones(getnbnode(fullModel),1)) ;
    else
        source = smooth(pb.model,source) ;
        qpCoord = getDiscon2Con(pb.model)*getDomainCoord(pb.model) ;
        source = relativeSort(source,qpCoord,getcoord(getnode(fullModel)),...
            getfemobjectoptions('tolerancepoint')) ;
    end
    
end

end