classdef QPConductivityAssembler
    
    properties( SetAccess = public, GetAccess = public )
        conductivity
        model % QPModel
        patterns
        patternsTable
        fields
        distributor
        distribution
        conductivityBounds
    end
    
    methods( Access = public )
        %% Constructor methods
        function assembler = QPConductivityAssembler(varargin)
            p = ImprovedInputParser();
            addParameter(p,'model',[]); % required
            addParameter(p,'patterns',[]);
            addParameter(p,'patternsTable',[]);
            addParameter(p,'fields',[]);
            addParameter(p,'distributor',[]);
            addParameter(p,'distribution',[]);
            addParameter(p,'conductivity',[]);
            parse(p,varargin{:}) ;
            assembler = passMatchedArgsToProperties(p,assembler);
            if isempty(getDistributor(assembler))
                proba = getcharin('probability',varargin,[0.9 0.1]) ;
                defaultDistributor = @(n) dealMultinomial(proba,n);
                assembler = setDistributor(assembler,defaultDistributor) ;
            end
            if isempty(getDistribution(assembler))
                assembler = setDistribution(assembler,...
                    distribute(assembler)) ;
            end
            pattern = getPatterns(assembler) ;
            if ~isempty(pattern) && isempty(getPatternsTable(assembler))
                table = eye(numel(pattern)) ;
                assembler = setPatternsTable(assembler,table) ;
            end
        end
        
        %% Get methods
        
        % Primary
        
        function conductivity = getConductivity(Assembler)
            conductivity = Assembler.conductivity ;
        end
        
        function model = getModel(Assembler)
            model = Assembler.model ;
        end
        
        function patterns = getPatterns(Assembler)
            patterns = Assembler.patterns ;
        end
        
        function patternsTable = getPatternsTable(Assembler)
            patternsTable = Assembler.patternsTable ;
        end
        
        function fields = getFields(Assembler)
            fields = Assembler.fields ;
        end
        
        function distributor = getDistributor(Assembler)
            distributor = Assembler.distributor ;
        end
        
        function distribution = getDistribution(Assembler)
            distribution = Assembler.distribution ;
        end
        
        function conductivityBounds = getConductivityBounds(Assembler,choice)
            if nargin==1 ; choice = 'all' ; end
            conductivityBounds = Assembler.conductivityBounds ;
            switch choice
                case {'all',1:2} % do nothing
                case {'min',1}
                    conductivityBounds = conductivityBounds(:,1) ;
                case {'max',2}
                    conductivityBounds = conductivityBounds(:,2) ;
                otherwise
                    error('choice is either "all", "min" or "max".')
            end
        end
        
        % Secondary
        
        function phasesNb = getPhasesNb(assembler)
            K = getConductivity(assembler) ;
            if ~isempty(K)
                phasesNb = K.space.dim(end) ; % actual pattern number
            else
                dist = getDistribution(assembler) ;
                phasesNb = sum(~cellfun(@isempty,dist)) ; % actual pattern number
                if ~phasesNb % distribution is empty
                    phasesNb = size(getFields(assembler),2); % max pattern number
                    if ~phasesNb % fields are empty
                        phasesNb = numel(distribute(assembler)) ;
                        % pattern number for a random realization
                    end
                end
            end
        end
        
        % Reference to model
        
        function order = getOrder(assembler)
            order = getOrder(getModel(assembler)) ;
        end
        
        function cellModel = getCellModel(Assembler,varargin)
            cellModel = getCellModel(getModel(Assembler),varargin{:}) ;
        end
        
        function cellNum = getCellNum(Assembler)
            cellNum = getCellNum(getModel(Assembler)) ;
        end
        
        function cellNb = getCellNb(Assembler)
            cellNb = getCellNb(getModel(Assembler)) ;
        end
        
        function cellSize = getCellSize(Assembler)
            cellSize = getCellSize(getModel(Assembler)) ;
        end
        
        function cellCoord = getCellCoord(Assembler)
            cellCoord = getCellCoord(getModel(Assembler)) ;
        end
        
        function tolSVD = getTolSVD(Assembler)
            tolSVD = getTolSVD(getModel(Assembler)) ;
        end
        
        function verbose = getVerbose(Assembler)
            verbose = getVerbose(getModel(Assembler)) ;
        end
        
        %% Set methods
        
        function Assembler = setConductivity(Assembler,newConductivity)
            Assembler.conductivity = newConductivity ;
        end
        
        function Assembler = setModel(Assembler,newModel)
            Assembler.model = newModel ;
        end
        
        function Assembler = setPatterns(Assembler,newPatterns,phaseNum)
            if nargin < 3
                phaseNum = 1:getPhasesNb(Assembler) ;
            end
            currentPatterns = getPatterns(Assembler) ;
            currentPatterns(phaseNum) = newPatterns ;
            Assembler.patterns = currentPatterns ;
        end
        
        function Assembler = setPatternsTable(Assembler,newPatternsTable)
            Assembler.patternsTable = newPatternsTable ;
        end
        
        function Assembler = setFields(Assembler,newFields,phaseNum)
            if nargin < 3
                phaseNum = [] ;
            end
            currentFields = getFields(Assembler) ;
            if isempty(currentFields) || isempty(phaseNum) % safety
                currentFields = newFields ;
            else
                currentFields(:,phaseNum) = newFields ;
            end
            Assembler.fields = currentFields ;
        end
        
        function Assembler = setDistributor(Assembler,newDistributor)
            Assembler.distributor = newDistributor ;
        end
        
        function Assembler = setDistribution(Assembler,newDistribution)
            Assembler.distribution = newDistribution ;
        end
        
        function Assembler = setConductivityBounds(Assembler,newBounds)
            if nargin==1
                newBounds = calcConductivityBounds(Assembler) ;
            end
            Assembler.conductivityBounds = newBounds ;
        end
        
        function Assembler = setTolSVD(Assembler,newTolSVD)
            newModel = setTolSVD(getModel(Assembler),newTolSVD) ;
            Assembler = setModel(Assembler,newModel) ;
        end
        
        function Assembler = setVerbose(Assembler,newVerbose)
            newModel = setVerbose(getModel(Assembler),newVerbose) ;
            Assembler = setModel(Assembler,newModel) ;
        end
        
        %% External methods signatures
        
        [Assembler,assemblerTime] = assemble(Assembler)
        
        assembler = assemblePatterns(assembler,patternsTable,patterns)
        
        conductivityBounds = calcConductivityBounds(assembler,method)
        
        distribution = distribute(Assembler)
        
        [] = ifprint(assembler,string)
        
        [Assembler,assemblerTime] = updateDistribution(Assembler,newDistribution)
        
        assembler = updateModel(assembler,model)
        
    end
    
    methods (Static)
        
        assembler = createRandom(varargin)
        
        assembler = testPatterns(varargin)
        
        assembler = fromConductivity(conductivity,model)
        
        assembler = homogeneous(value,model)
        
    end
end