function [operator,time] = assembleStabilisationOperator(assembler)
% [operator,time] = assembleStabilisationOperator(assembler)

assemblerClock = tic ;

% Get elementary operators
microOperators = getMicroOperators(assembler) ;
microStabilisationOperators = microOperators.stabilisation;
mesoOperators = getMesoOperators(assembler) ;
mesoStabilisationOperators = mesoOperators.stabilisation;
penalty = getPenalty(assembler) ;

% Build Tspace and core
% Each dimension's operators stored as single column cell array:
mesoStabilisationOperators = cellfun(@(x) reshape(x,[],1),...
    mesoStabilisationOperators,'UniformOutput',false) ;
space = TSpaceOperators([mesoStabilisationOperators(:) ; ...
    {microStabilisationOperators(:)}]) ;
switch getOrder(assembler)
    case 2
        core = zeros(8,8) ;
        core(1,1) = 1 ; % l({}^t\ki^1\odot{}^t\omega) \otimes M_0^e1
        core(2,2) = 1 ; % l(\ki^1\odot\omega) \otimes M_0^-e1
        core(3,3) = -1 ; % -{}^t\ki^1\odot{}^t\omega) \otimes M_1^e1
        core(4,4) = -1 ; % -\ki^1\odot\omega) \otimes M_1^-e1
        core(5,5) = 1 ; % l({}^t\ki^2\odot{}^t\omega) \otimes M_0^e2
        core(6,6) = 1 ; % l(\ki^2\odot\omega) \otimes M_0^-e2
        core(7,7) = -1 ; % -{}^t\ki^2\odot{}^t\omega) \otimes M_1^e2
        core(8,8) = -1 ; % -\ki^2\odot\omega) \otimes M_1^-e2
        % Core could be defined as DiagonalTensor, but this is easier 
        % to edit and more uniform across both structures.
    case 3
        core = zeros(5,5,5) ;
        core(1,5,1) = 1 ;% l({}^t\ki^I\odot{}^t\omega^I_n)\otimes Id^J\odot\omega^J_n\otimes M_0^e1
        core(2,5,2) = 1 ;% l(\ki^I\odot\omega^I_n)\otimes Id^J\odot\omega^J_n\otimes M_0^-e1
        core(3,5,3) = -1 ;% -{}^t\ki^I\odot{}^t\omega^I_n\otimes Id^J\odot\omega^J_n\otimes M_1^e1
        core(4,5,4) = -1 ;% -\ki^I\odot\omega^I_n\otimes Id^J\odot\omega^J_n\otimes M_1^-e1
        core(5,1,5) = 1 ;% Id^I\odot\omega^I_n\otimes l({}^t\ki^J\odot{}^t\omega^J_n)\otimes M_0^e2
        core(5,2,6) = 1 ;% Id^I\odot\omega^I_n\otimes l(\ki^J\odot\omega^J_n)\otimes M_0^-e2
        core(5,3,7) = -1 ;% -Id^I\odot\omega^I_n\otimes{}^t\ki^J\odot{}^t\omega^J_n\otimes M_1^e2
        core(5,4,8) = -1 ;% -Id^I\odot\omega^I_n\otimes\ki^J\odot\omega^J_n\otimes M_1^-e2
end
%TODO: exploit core sparsity ?
weightRank = space.dim./size(core) ;% weight operator's rank
weightRank = weightRank(1:end-1) ;% no rank along order associated to Y
weightCore = eye(weightRank) ;% canonical format
core = superkron(weightCore,core) ;% Operators in TSpace grouped by same \omega_n
core = FullTensor(penalty*core) ;

% Assemble into tensor
operator = TuckerLikeTensor(core,space) ;

time = toc(assemblerClock) ;
end