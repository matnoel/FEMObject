function [operator,time] = assembleDiffusionOperator(assembler)
% [operator,time] = assembleDiffusionOperator(assembler)

assemblerClock = tic ;

% Get elementary operators
microOperators = getMicroOperators(assembler) ;
microDiffusionOperators = microOperators.diffusion ;
mesoOperators = getMesoOperators(assembler) ;
mesoDiffusionOperators = mesoOperators.diffusion ;
conductivity = getConductivity(assembler) ;

% Build Tspace and core
% Each dimension's operators stored as single column cell array:
mesoDiffusionOperators = cellfun(@(x) reshape(x,[],1),...
    mesoDiffusionOperators,'UniformOutput',false) ;
space = TSpaceOperators([mesoDiffusionOperators(:) ; ...
    {microDiffusionOperators(:)}]) ;
core = conductivity.core ;

% Assemble into tensor
operator = TuckerLikeTensor(core,space) ;

time = toc(assemblerClock) ;

end