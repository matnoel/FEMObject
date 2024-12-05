classdef QPDiffusionAssembler
    
    properties( SetAccess = public, GetAccess = public )
        conductivityAssembler
        bc
        source
        patches
        constantNullification
        penalty
        microOperators
        mesoOperators
        lhsOperator
        rhsOperator
        stabilisationWeights
        averageWeights
        useStabilisationWeights
        useAverageWeights
        useCompression
    end
    
    methods( Access = public )
        %% Constructor methods
        function assembler = QPDiffusionAssembler(varargin)
            p = ImprovedInputParser();
            addParameter(p,'conductivityAssembler',[]);
            addParameter(p,'bc','PBC');
            addParameter(p,'source',[]) ;
            addParameter(p,'patches',{});
            addParameter(p,'penalty',[]);
            addParameter(p,'constantNullification','full');
            addParameter(p,'useStabilisationWeights',false);
            addParameter(p,'useAverageWeights',false);
            addParameter(p,'useCompression',false);
            addParameter(p,'microOperators',[]);
            addParameter(p,'mesoOperators',[]);
            addParameter(p,'stabilisationWeights',[]);
            addParameter(p,'averageWeights',[]);
            parse(p,varargin{:}) ;
            assembler = passMatchedArgsToProperties(p,assembler);
            % Initialize
            if isempty(getSource(assembler))
                assembler = setSource(assembler,0) ;
            else
                assembler = setSource(assembler) ;
            end
            % Format bc attribute as cell array of QPBC instances
            rawBC = getBC(assembler) ;
            if ischar(rawBC) && strcmp('PBC',rawBC)
                assembler = setBC(assembler,...
                    {QPBC('type',5,'model',getModel(assembler))}) ;
            elseif isnumeric(rawBC)
                assembler = setBC(assembler,...
                    {QPBC('type',1,'model',getModel(assembler),'value',rawBC)}) ;
            elseif isa(rawBC,'QPBC')
                assembler = setBC(assembler,{rawBC}) ;
            end
            % Compute weights and penalty
            if ~isempty(getConductivityBounds(assembler)) % safety
                if isempty(getPenalty(assembler))
                    assembler = setPenalty(assembler) ;
                end
                if isempty(getStabilisationWeights(assembler)) &&...
                    getUseStabilisationWeights(assembler)
                    assembler = setStabilisationWeights(assembler) ;
                end
                if isempty(getAverageWeights(assembler)) &&...
                        getUseAverageWeights(assembler)
                    assembler = setAverageWeights(assembler) ;
                end
            end
        end
        
        %% "Get" methods
        
        % Primary
        
        function conductivityAssembler = getConductivityAssembler(assembler)
            conductivityAssembler = assembler.conductivityAssembler ;
        end
        
        function bc = getBC(assembler)
            bc = assembler.bc ;
        end
        
        function source = getSource(assembler)
            source = assembler.source ;
        end
        
        function patches = getPatches(assembler)
            patches = assembler.patches ;
        end
        
        function penalty = getPenalty(assembler)
            penalty = assembler.penalty ;
        end
        
        function constantNullification = getConstantNullification(assembler)
            constantNullification = assembler.constantNullification ;
        end
        
        function microOperators = getMicroOperators(assembler)
            microOperators = assembler.microOperators ;
        end
        
        function mesoOperators = getMesoOperators(assembler)
            mesoOperators = assembler.mesoOperators ;
        end
        
        function sWeights = getStabilisationWeights(assembler,order)
            if nargin == 1
                order = 1:getOrder(assembler)-1 ;
            end
            sWeights = assembler.stabilisationWeights ;
            if ~isempty(sWeights) % safety
                sWeights = sWeights{order} ;
            end
        end
        
        function aWeights = getAverageWeights(assembler,order)
            if nargin == 1
                order = 1:getOrder(assembler)-1 ;
            end
            aWeights = assembler.averageWeights ;
            if ~isempty(aWeights) % safety
                aWeights = aWeights{order} ;
            end
        end
        
        function lhsOperator = getLHSOperator(assembler)
            lhsOperator = assembler.lhsOperator ;
        end
        
        function rhsOperator = getRHSOperator(assembler)
            rhsOperator = assembler.rhsOperator ;
        end
        
        function useAverageWeights = getUseAverageWeights(assembler)
            useAverageWeights = assembler.useAverageWeights ;
        end
        
        function useStabilisationWeights = getUseStabilisationWeights(assembler)
            useStabilisationWeights = assembler.useStabilisationWeights ;
        end
        
        function useCompression = getUseCompression(assembler)
            useCompression = assembler.useCompression ;
        end
        
        % Secondary
        
        function patternNb = getConductivityPatternNb(assembler)
            patternNb = getPhasesNb(getConductivityAssembler(assembler)) ;
        end
        
        % References to conductivity assembler
        
        function conductivity = getConductivity(assembler)
            conductivity = getConductivity(getConductivityAssembler(assembler)) ;
        end
        
        function conductivityBounds = getConductivityBounds(assembler,choice)
            if nargin < 2
                choice = 'all' ;
            end
            conductivityBounds = getConductivityBounds(...
                getConductivityAssembler(assembler),choice) ;
        end
        
        function distribution = getDistribution(assembler)
            distribution = getDistribution(getConductivityAssembler(assembler)) ;
        end
        
        % References to model
        
        function order = getOrder(assembler)
            order = getOrder(getModel(assembler)) ;
        end
        
        function model = getModel(assembler)
            model = getModel(getConductivityAssembler(assembler)) ;
        end
        
        function traceConstantL2 = getTraceConstantL2(assembler)
            traceConstantL2 = getTraceConstantL2(getModel(assembler)) ;
        end
        
        function tolSVD = getTolSVD(assembler)
            tolSVD = getTolSVD(getModel(assembler)) ;
        end
        
        function verbose = getVerbose(assembler)
            verbose = getVerbose(getModel(assembler)) ;
        end
        
        function cellModel = getCellModel(assembler,varargin)
            cellModel = getCellModel(getModel(assembler),varargin{:}) ;
        end
        
        function cellNum = getCellNum(assembler)
            cellNum = getCellNum(getModel(assembler)) ;
        end
        
        function cellNb = getCellNb(assembler)
            cellNb = getCellNb(getModel(assembler)) ;
        end
        
        function cellSize = getCellSize(assembler)
            cellSize = getCellSize(getModel(assembler)) ;
        end
        
        function cellDomain = getCellDomain(assembler)
            cellDomain = getCellDomain(getModel(assembler)) ;
        end
        
        function coord = getCoord(assembler,varargin)
            coord = getCoord(getModel(assembler),varargin{:}) ;
        end
        
        function nbCellDoF = getNbCellDoF(assembler)
            nbCellDoF = getNbCellDoF(getModel(assembler));
        end
        
        function nbDomainDoF = getNbDomainDoF(assembler)
            nbDomainDoF = getNbDomainDoF(getModel(assembler)) ;
        end
        
        % References to patches
        
        function patchesCells = getPatchesCells(assembler)
            ps = getPatches(assembler) ;
            patchesCells = [] ;
            if ~isempty(ps)
                patchesCells = getCellList(ps) ;
            end
        end
        
        function patchesNb = getPatchesNb(assembler)
            patchesNb = numel(getPatches(assembler)) ;
        end
        
        %% "Set" methods
        
        function assembler = setConductivityAssembler(assembler,conductivityAssembler)
            assembler.conductivityAssembler = conductivityAssembler ;
        end
        
        function assembler = setPenalty(assembler,newPenalty)
            if nargin < 2
                newPenalty = calcPenalty(assembler) ;
            end
                assembler.penalty = newPenalty ;
                assembler = setBCPenalty(assembler,newPenalty) ;
        end
        
        function assembler = setBC(assembler,newBC)
            if ~iscell(newBC)
                newBC = {newBC} ;
            end
                assembler.bc = newBC ;
        end
                
        function assembler = setConstantNullification(assembler,newCstNull)
            assembler.constantNullification = newCstNull ;
        end
        
        function assembler = setSource(assembler,newSource)
            if nargin == 1
                newSource = getSource(assembler) ;
            end
            assembler.source = formatSource(assembler,newSource) ;
        end
        
        function assembler = setPatches(assembler,patches)
            assembler.patches = patches ;
        end
        
        function assembler = setUseStabilisationWeights(assembler,newValue)
            assembler.useStabilisationWeights = newValue ;
        end
        
        function assembler = setUseAverageWeights(assembler,newValue)
            assembler.useAverageWeights = newValue ;
        end
        
        function assembler = setStabilisationWeights(assembler,newStabilisationWeights)
            if nargin == 1
                newStabilisationWeights = calcStabilisationWeights(assembler) ;
            end
            assembler.stabilisationWeights = newStabilisationWeights ;
        end
        
        function assembler = setAverageWeights(assembler,newAverageWeights)
            if nargin == 1
                newAverageWeights = calcAverageWeights(assembler) ;
            end
            assembler.averageWeights = newAverageWeights ;
        end
        
        function assembler = setMesoOperators(assembler,mesoOperators)
            assembler.mesoOperators = mesoOperators ;
        end
        
        function assembler = setMicroOperators(assembler,microOperators)
            assembler.microOperators = microOperators ;
        end
        
        function assembler = setLHSOperator(assembler,lhsOperator)
            assembler.lhsOperator = lhsOperator ;
        end
        
        function assembler = setRHSOperator(assembler,rhsOperator)
            assembler.rhsOperator = rhsOperator ;
        end
        
        function assembler = setUseCompression(assembler,newUseCompression)
            assembler.useCompression = newUseCompression ;
        end
        
        % Secondary
        
        function assembler = addBC(assembler,newBC)
            if ~iscell(newBC)
                newBC = {newBC} ;
            end
            bc = getBC(assembler) ;
            assembler = setBC(assembler,[bc(:) ; newBC]) ;
        end
                
        % References
        
        function assembler = setVerbose(assembler,newVerbose)
            cA = setVerbose(getConductivityAssembler(assembler),newVerbose) ;
            assembler = setConductivityAssembler(assembler,cA) ;
        end
        
        function assembler = setBCPenalty(assembler,newPenalty)
            if nargin < 2
                newPenalty = getPenalty(assembler) ;
            end
            bc = getBC(assembler) ;
            for i = 1:numel(bc)
                if ismember(getType(bc{i}),[1 4]) % Dirichlet or Cauchy
                    bc{i} = setPenalty(bc{i},newPenalty) ;
                end
            end
            assembler = setBC(assembler,bc) ;
        end
        
        function assembler = addPatch(assembler,patch)
            patches = getPatches(assembler) ;
            patches = addPatch(patches,patch) ;
            assembler = setPatches(assembler,patches) ;
        end
        
        %% External methods signatures
        
        % Assembling methods
        
        [assembler,assemblerTime] = assemble(assembler,orders2Assemble)
        
        [lhs,rhs,time] = assembleBCOperator(assembler)
        
        [operator,time] = assembleConsistencyOperator(assembler)
        
        [operator,time] = assembleCstNullOperator(assembler)
        
        [operator,time] = assembleDiffusionOperator(assembler)
        
        [assembler,assemblerTime] = assembleMesoOperators(assembler)
        
        mesoOperators = assembleMesoOperatorsDirection(assembler,dir)
        
        [assembler,assemblerTime] = assembleMicroOperators(assembler)
        
        [operator,time] = assembleRHSOperator(assembler)
        
        [operator,time] = assembleStabilisationOperator(assembler)
        
        % Miscellaneous
        
        weights = calcAverageWeights(assembler)
        
        flux = calcBoundaryFlux(assembler,v,bCells)
        
        cellConn = calcMesoConnectivity(assembler,direction)
        
        penalty = calcPenalty(assembler)
        
        weights = calcStabilisationWeights(assembler)
        
        [assembler, compressionTime] = compressOperators(assembler,tolerance)
        
        source = formatSource(assembler,source)
        
        [] = ifprint(assembler,string)
        
        bool = isPeriodic(assembler,direction)
        
        [assembler,time] = updateBC(assembler,bc,bcNum)
        
        diffusionAssembler = updateConductivityAssembler(diffusionAssembler,conductivityAssembler)
        
        assembler = updateModel(assembler,model)
        
        [assembler,time] = updateSource(assembler,source)
        
    end
    
    methods (Static)
        
        assembler = createRandom(varargin)
        
    end
end