function source = formatSource(assembler,source)
% source = formatSource(assembler,source)

order = getOrder(assembler) ;

switch class(source)
    
    case 'cell'
        switch class(source{1})
            case 'double'
                space = TSpaceVectors(source(:)) ;
            case 'cell'
                space = TSpaceOperators(source(:)) ;
                warning('QPDiffusionAssembler.formatSource: source should not be a TSpaceOperators')
            otherwise
                error('QPDiffusionAssembler.formatSource: unknown source class in cell')
        end
        core = DiagonalTensor(ones(min(space.dim),1),order) ;
        source = TuckerLikeTensor(core,space) ;
        source = formatSource(assembler,source) ;
        
    case 'double'
        nbCellDoF = getNbCellDoF(assembler);
        nbDomainDoF = getNbDomainDoF(assembler) ;
        nbTotalDoF = getNbTotalDoF(getModel(assembler)) ;
        switch numel(source)
            case 0 % source is empty
                source = 0 ;
                source = formatSource(assembler,source) ;
            case 1 % source is constant
                source = source*ones(nbCellDoF,1) ;
                source = formatSource(assembler,source) ;
            case nbCellDoF % source is scalar field over cell
                model = getModel(assembler) ;
                cellNb = getCellNb(assembler) ;
                mesoIndic = mesoIndicator(model,(1:cellNb)',0) ;
                source = [mesoIndic(:) ; {source}] ;
                source = formatSource(assembler,source) ;
            case {nbDomainDoF,nbDomainDoF^2,nbTotalDoF} % source is scalar field or operator over whole domain
                source = tensorize(getModel(assembler),source) ;
                source = formatSource(assembler,source) ;
            case nbCellDoF^2 % source is operator over cell (proper shape assumed)
                warning('Source should not be operator')
                cellNb = getCellNb(assembler) ;
                if order==2
                    indic = coord2indicator(cellNb,1:cellNb) ;
                elseif order==3
                    indic = coord2indicator(getCellNum(assembler),1:cellNb) ;
                end
                % format like TSpaceOperators
                indic = cellfun(@(x) {diag(x)},indic(:),'UniformOutput',false) ;
                source = [indic ; {{source}}] ;
                source = formatSource(assembler,source) ;
            otherwise
                error('Source size mismatch.')
        end
        
    case {'FullTensor','DiagonalTensor'}
        source = formatSource(assembler,double(source)) ;
        
    case 'CanonicalTensor'
        source = TuckerLikeTensor(source.core,source.space) ;
        source = formatSource(assembler,source) ;
        
    case 'TuckerLikeTensor'
        if source.order == order
            return
            % End of recursion: intended format.
        else
            source = formatSource(assembler,untensorize(...
                getModel(assembler),source)) ;
            %TODO: brutal (and costly)
        end
        
    case {'TSpaceVectors','TSpaceOperators'}
        source = formatSource(assembler,source.spaces) ;
        
    case 'function_handle' % function to apply to relative coordinates (costly)
        coord = getDomainCoord(getModel(assembler)) ;
        cellSize = getCellSize(assembler) ;
        source = source(coord*diag(cellSize.^-1)) ; % apply to *relative* coordinates
        source = formatSource(assembler,source) ;
        
    case 'char'
        assert(strcmp('corrector1',source) ...
               || strcmp('corrector2',source) ...
               || strcmp('correctors',source), ... % FE solver only
            'Unknown source string') ;
        % End of recursion: will be processed in assembleRHSOperator
        
    otherwise
        error('QPDiffusionAssembler.formatSource: unknown source class')
end

end