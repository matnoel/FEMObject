function [operator,time] = assembleConsistencyOperator(assembler)
% [operator,time] = assembleConsistencyOperator(assembler)

assemblerClock = tic ;

% Get elementary operators
microOperators = getMicroOperators(assembler) ;
microConsistencyOperators = microOperators.consistency;
mesoOperators = getMesoOperators(assembler) ;
mesoConsistencyOperators = mesoOperators.consistency;
conductivity = getConductivity(assembler) ;

% Build Tspace and core
% Each dimension's operators stored as single column cell array:
mesoConsistencyOperators = cellfun(@(x) reshape(x,[],1),...
    mesoConsistencyOperators,'UniformOutput',false) ;
space = TSpaceOperators([mesoConsistencyOperators(:) ; ...
    {microConsistencyOperators(:)}]) ;
switch getOrder(assembler)
    case 2
        core = zeros(8,8) ;
        core(1,1) = 1 ; % l({}^t\ki^1\odot K^I_n) \otimes N_0^e1[K^Y_n]
        core(2,2) = 1 ; % l(\ki^1\odot K^I_n) \otimes N_0^-e1[K^Y_n]
        core(3,3) = -1 ; % -{}^t\ki^1\odot K^I_n \otimes N_1^e1[K^Y_n]
        core(4,4) = -1 ; % -\ki^1\odot K^I_n) \otimes N_1^-e1[K^Y_n]
        core(5,5) = 1 ; % l({}^t\ki^2\odot K^I_n) \otimes N_0^e2[K^Y_n]
        core(6,6) = 1 ; % l(\ki^2\odot K^I_n) \otimes N_0^-e2[K^Y_n]
        core(7,7) = -1 ; % -{}^t\ki^2\odot K^I_n \otimes N_1^e2[K^Y_n]
        core(8,8) = -1 ; % -\ki^2\odot K^I_n) \otimes N_1^-e2[K^Y_n]
        % Core could be defined as DiagonalTensor, but this is easier 
        % to edit and more uniform across both structures.
    case 3
        core = zeros(5,5,5) ;
        core(1,5,1) = 1 ;% l({}^t\ki^I\odot{}^t K^I_n)\otimes Id^J\odot K^J_n\otimes N_0^e1
        core(2,5,2) = 1 ;% l(\ki^I\odot K^I_n)\otimes Id^J\odot K^J_n\otimes N_0^-e1
        core(3,5,3) = -1 ;% -{}^t\ki^I\odot{}^t K^I_n\otimes Id^J\odot K^J_n\otimes N_1^e1
        core(4,5,4) = -1 ;% -\ki^I\odot K^I_n\otimes Id^J\odot K^J_n\otimes N_1^-e1
        core(5,1,5) = 1 ;% Id^I\odot K^I_n\otimes l({}^t\ki^J\odot{}^t K^J_n)\otimes N_0^e2
        core(5,2,6) = 1 ;% Id^I\odot K^I_n\otimes l(\ki^J\odot K^J_n)\otimes N_0^-e2
        core(5,3,7) = -1 ;% -Id^I\odot K^I_n\otimes{}^t\ki^J\odot{}^t K^J_n\otimes N_1^e2
        core(5,4,8) = -1 ;% -Id^I\odot K^I_n\otimes\ki^J\odot K^J_n\otimes N_1^-e2
end
%TODO: exploit core sparsity ?
core = superkron(double(conductivity.core),core) ;
weightRank = space.dim./size(core) ;% weight operator's rank
weightRank = weightRank(1:end-1) ;% no rank along order associated to Y
weightCore = eye(weightRank) ;% canonical format
core = superkron(weightCore,core) ;% Operators in TSpace grouped by same \omega_n
core = FullTensor(core) ;

% Assemble into tensor
operator = TuckerLikeTensor(core,space) ;

time = toc(assemblerClock) ;

end