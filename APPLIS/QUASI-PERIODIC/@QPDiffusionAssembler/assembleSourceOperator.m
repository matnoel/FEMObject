function [sourceOperator,time] = assembleSourceOperator(assembler)
% [operator,time] = assembleSourceOperator(assembler)

assemblerClock = tic ;

source = getSource(assembler) ;

switch class(source)
    case 'char' % corrector source term
        
        % Get operators
        diffusionOperator = assembleDiffusionOperator(assembler) ;
        consistencyOperator = assembleConsistencyOperator(assembler) ;
        
        % Get coordinates
        excepted = getPatchesCells(assembler) ;
        x = getCoord(assembler,'except',excepted) ;
        if strcmp('corrector1',source)
            x = x{1} ;
        elseif strcmp('corrector2',source)
            x = x{2} ;
        else
            error('Unknown source string')
        end
        
        % Assemble and store RHS operator
        sourceOperator = diffusionOperator - consistencyOperator ;
        sourceOperator = sourceOperator*x ;
        
    case {'TuckerLikeTensor','FullTensor','AlgebraicTensor','TTTensor','CanonicalTensor'}
        
        if isa(source.space,'TSpaceVectors') % standard assembling
            % Get elementary operators
            microOperators = getMicroOperators(assembler) ;
            microRHSOperators = microOperators.rhs ;
            mesoOperators = getMesoOperators(assembler) ;
            mesoRHSOperators = mesoOperators.rhs ;
            
            % Build Tspace and core
            space = TSpaceVectors([mesoRHSOperators(:) ; ...
                {microRHSOperators}]) ;
            core = source.core ;
            
            % Assemble into tensor
            sourceOperator = TuckerLikeTensor(core,space) ;
            
        elseif isa(source.space,'TSpaceOperators') % copy source
            sourceOperator = source ;
            warning('Source term copied as RHS operator')
        end
        
        if ~isempty(getPatchesCells(assembler))
            nullSource = restrictTensor(getModel(assembler),sourceOperator,...
                getPatchesCells(assembler),0) ; % 0 : no orthogonalization
            sourceOperator = sourceOperator - nullSource ;
        end
        
    otherwise
        error('Unknown source class')
end

time = toc(assemblerClock) ;
end