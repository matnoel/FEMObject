function [operator,time] = assembleCstNullOperator(assembler)
% [operator,time] = assembleCstNullOperator(assembler)

assemblerClock = tic ;

% If not all BC are either Neumann or Periodic, return 0
if any(~ismember(cellfun(@getType,getBC(assembler)),[2 5]))
    tsz = tensorSize(getModel(assembler)) ;
    operator = TuckerLikeTensor.zeros(tsz,tsz) ;
    time = toc(assemblerClock) ;
    return
end

% Get elementary operators
microOperators = getMicroOperators(assembler) ;
microCstNullOperator = microOperators.cstNull ;
mesoOperators = getMesoOperators(assembler) ;
mesoCstNullOperator = mesoOperators.cstNull ;

% Build Tspace and core
% Each dimension's operators stored as single column cell array:
mesoCstNullOperator = cellfun(@(x) reshape(x,[],1),...
    mesoCstNullOperator,'UniformOutput',false) ;
space = TSpaceOperators([mesoCstNullOperator(:) ; ...
    {microCstNullOperator(:)}]) ;
core = DiagonalTensor(ones(min(space.dim),1),getOrder(assembler)) ;

% Assemble into tensor
operator = TuckerLikeTensor(core,space) ;

time = toc(assemblerClock) ;
end