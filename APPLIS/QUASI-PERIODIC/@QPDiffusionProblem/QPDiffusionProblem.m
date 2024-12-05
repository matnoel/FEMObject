classdef QPDiffusionProblem
    
    properties
        operatorAssembler
        greedySolver
        updater
        tolerance
    end
    
    methods( Access = public )
        %% Constructor methods
        function pb = QPDiffusionProblem(varargin)
            p = ImprovedInputParser();
            addParameter(p,'operatorAssembler',[]);
            addParameter(p,'tolerance',1e-3);
            addParameter(p,'greedySolver',[]);
            addParameter(p,'updater',[]);
            parse(p,varargin{:});
            pb = passMatchedArgsToProperties(p,pb);
            if ~isempty(getOperatorAssembler(pb)) && ~isempty(getLHSOperator(pb))
                pb = setDefaults(pb,false) ;
            end
        end
        
        function pb = setDefaults(pb,overwrite)
            if nargin == 1
                overwrite = true ;
            end
            if overwrite || isempty(getUpdater(pb))
                pb = setUpdater(pb,defaultUpdater(pb)) ;
            end
            if overwrite || isempty(getGreedySolver(pb))
                pb = setGreedySolver(pb,defaultGreedySolver(pb)) ;
            end
        end
        
        function solver = defaultLocalSolver(pb)
            solver = RankOneALSLinearSolver(getLHSOperator(pb),...
                getRHSOperator(pb),...
                'maxIterations',10,...
                'stagnation',getTolerance(pb)/100,...
                'display',false,...
                'x0',[],...
                'tolerance',1e-3) ; % tolerance is unused
        end
        
        function updater = defaultUpdater(pb)
            tol = getTolerance(pb) ;
            updater = TuckerLikeTensorALSLinearSolver(getLHSOperator(pb),...
                getRHSOperator(pb),...
                'tolerance',tol,...
                'maxIterationsTSpace',1,...
                'stagnationTSpace',min(1e-6,tol/100),...
                'display',getVerbose(pb)) ;
        end
        
        function greedy = defaultGreedySolver(pb)
            tol = getTolerance(pb) ;
            up = getUpdater(pb) ;
            greedy = GreedyLinearSolver(getLHSOperator(pb),...
                getRHSOperator(pb),...
                'localSolver',defaultLocalSolver(pb),...
                'update',@(x)updateTSpaceByALS(up,x),...
                'x0',[],...
                'tolerance',tol, ...
                'maxIterations',getCellNb(pb),...
                'stagnation',min(1e-6,tol/100),...
                'minimizeResidual',false,... % ?
                'checkResidual',1,... % Set to N to check residual every N steps
                'display',getVerbose(pb));
        end
                
        %% "Get" methods
        
        % Primary
        
        function greedySolver = getGreedySolver(pb)
            greedySolver = pb.greedySolver ;
        end
        
        function operatorAssembler = getOperatorAssembler(pb)
            operatorAssembler = pb.operatorAssembler ;
        end
        
        function tolerance = getTolerance(pb)
            tolerance = pb.tolerance ;
        end
        
        % Secondary
                
        function localSolver = getLocalSolver(pb)
            greedy = getGreedySolver(pb) ;
            if isempty(greedy)
                localSolver = [] ;
            else
                localSolver = greedy.localSolver ;
            end
        end
        
        function updater = getUpdater(pb)
            updater = pb.updater ;
        end
        
        function initialPoint = getInitialPoint(pb)
            greedy = getGreedySolver(pb) ;
            if isempty(greedy)
                initialPoint = [] ;
            else
                initialPoint = greedy.x0 ;
            end
        end
        
        % References to operator assembler
        
        function model = getModel(pb)
            model = getModel(getOperatorAssembler(pb)) ;
        end
        
        function lhsOperator = getLHSOperator(pb)
            lhsOperator = getLHSOperator(getOperatorAssembler(pb)) ;
        end
        
        function rhsOperator = getRHSOperator(pb)
            rhsOperator = getRHSOperator(getOperatorAssembler(pb)) ;
        end
        
        function microOp = getMicroOperators(pb)
            microOp = getMicroOperators(getOperatorAssembler(pb)) ;
        end
        
        function verbose = getVerbose(pb)
            verbose = getVerbose(getOperatorAssembler(pb)) ;
        end
        
        function conductivityAssembler = getConductivityAssembler(pb)
            conductivityAssembler = getConductivityAssembler(getOperatorAssembler(pb)) ;
        end
        
        function distribution = getDistribution(pb)
           distribution = getDistribution(getOperatorAssembler(pb)) ;
        end
        
        function conductivity = getConductivity(pb)
            conductivity = getConductivity(getOperatorAssembler(pb)) ;
        end
        
        function BC = getBC(pb)
            BC = getBC(getOperatorAssembler(pb)) ;
        end
        
        function source = getSource(pb)
            source = getSource(getOperatorAssembler(pb)) ;
        end
        
        function cstNull = getConstantNullification(pb)
            cstNull = getConstantNullification(getOperatorAssembler(pb)) ;
        end
        
        function conductivityBounds = getConductivityBounds(pb)
            conductivityBounds = getConductivityBounds(getOperatorAssembler(pb)) ;
        end
        
        function patches = getPatches(pb)
            patches = getPatches(getOperatorAssembler(pb)) ;
        end
        
        function tolSVD = getTolSVD(pb)
            tolSVD = getTolSVD(getOperatorAssembler(pb)) ;
        end
        
        % References to model
        
        function order = getOrder(pb)
            order = getOrder(getModel(pb)) ;
        end
        
        function cellNb = getCellNb(pb)
            cellNb = getCellNb(getModel(pb)) ;
        end
        
        function cellNum = getCellNum(pb)
            cellNum = getCellNum(getModel(pb)) ;
        end
        
        function cellSize = getCellSize(pb)
            cellSize = getCellSize(getModel(pb)) ;
        end
        
        function elementSize = getElementSize(pb)
            elementSize = getElemenSize(getModel(pb)) ;
        end
        
        function cellModel = getCellModel(pb,varargin)
            cellModel = getCellModel(getModel(pb),varargin{:}) ;
        end
        
        function cellDomain= getCellDomain(pb)
            cellDomain = getCellDomain(getModel(pb)) ;
        end
        
        function nbCellDoF = getNbCellDoF(pb) 
            nbCellDoF = getNbCellDoF(getModel(pb)) ;
        end
        
        function nbDomainDoF = getNbDomainDoF(pb) 
            nbDomainDoF = getNbDomainDoF(getModel(pb)) ;
        end
        
        % References to conductivity assembler
        
        function distributor = getDistributor(pb)
            distributor = getDistributor(getConductivityAssembler(pb)) ;
        end
        
        %% "Set" methods
        
        % Primary
        
        function pb = setGreedySolver(pb,newGreedySolver)
            if nargin == 1
                newGreedySolver = defaultGreedySolver(pb) ;
            end
            pb.greedySolver = newGreedySolver ;
        end
        
        function pb = setOperatorAssembler(pb,newOperatorAssembler)
            pb.operatorAssembler = newOperatorAssembler ;
        end
        
        function pb = setTolerance(pb,newTolerance)
            pb.tolerance = newTolerance ;
            greedy = getGreedySolver(pb) ;
            if ~isempty(greedy)
                greedy.tolerance = newTolerance ;
                pb = setGreedySolver(pb,greedy) ;
            end
        end
        
        function pb = setUpdater(pb,updater)
            if nargin == 1
                updater = defaultSolutionUpdater(pb) ;
            end
            pb.updater = updater ;
        end
        
        % References to operator assembler
        
        function pb = setBC(pb,newBC)
            opAss = getOperatorAssembler(pb) ;
            opAss = setBC(opAss,newBC) ;
            pb = setOperatorAssembler(pb,opAss) ;
        end
        
        function pb = setSource(pb,newSource)
            opAss = getOperatorAssembler(pb) ;
            opAss = setSource(opAss,newSource) ;
            pb = setOperatorAssembler(pb,opAss) ;
        end
                
        function pb = setVerbose(pb,newVerbose)
            oA = setVerbose(getOperatorAssembler(pb),newVerbose) ;
            pb = setOperatorAssembler(pb,oA) ;
        end
        
        % References to greedy solver
        
        function pb = setLocalSolver(pb,localSolver)
            if nargin == 1
                localSolver = defaultLocalSolver(pb) ;
            end
            greedy = getGreedySolver(pb) ;
            assert(~isempty(greedy),'No greedySolver') ;
            greedy.localSolver = localSolver ;
            pb = setGreedySolver(pb,greedy) ;
        end
        
        function pb = setInitialPoint(pb,initialPoint)
            greedy = getGreedySolver(pb) ;
            assert(~isempty(greedy),'No greedySolver') ;
            greedy.x0 = initialPoint ;
            greedy.initialUpdate = true ;
            pb = setGreedySolver(pb,greedy) ;
        end
        
        % References to conductivity assembler
        
        function pb = setDistributor(pb,newDistributor)
            condAss = getConductivityAssembler(getOperatorAssembler(pb)) ;
            condAss = setDistributor(condAss,newDistributor) ;
            pb = updateConductivityAssembler(pb,condAss) ;
        end
        
        function pb = setConductivity(pb,newConductivity)
        % pb = setConductivity(pb,newConductivity) 
        % WARNING: You should probably use updateConductivity instead.
            oA = getOperatorAssembler(pb) ;
            cA = getConductivityAssembler(oA) ;
            cA = setConductivity(cA,newConductivity) ;
            oA = setConductivityAssembler(oA,cA) ;
            pb = setOperatorAssembler(pb,oA) ; 
        end
        
        %% Display methods (model ?)
        
%         function [] = plotOutput(pb,output,name)
%         end
%         
%         function [] = plotRank(pb,output)
%         end
%         
%         function [] = plotTime(pb,output)
%         end
                
        %% External methods signatures
        
        correctedConductivity = applyCorrector(pb,corrector)
        
        [pb,assemblerTime] = assembleOperators(pb,orders2assemble)
        
        [errorVal,diff] = compare2FE(pb,tensQP,matFE,errorIndicator,varargin)
        
        [Khomo,output] = homogenise(pb,solverChoice,iterMax,stagnationTol,K)
        
        [] = ifprint(pb,string)
        
        bool = isPeriodic(pb,direction)
        
        [solutions,outputs,pb] = multiSolve(pb,nb)
        
        [solutions,outputs] = multiSolveFEM(pb,distributions,mesherMethod)
        
        [solution,output,pb] = solveDD(pb)
        
        [solution,output] = solveFEM(pb,useMultiplier,meshMethod)
        
        [solution,output,pb] = solve(pb,orders2assemble)
        
        [pb,time] = updateBC(pb,bc,bcNum)
        
        pb = updateConductivity(pb,K,distribution)
        
        pb = updateConductivityAssembler(pb,conductivityAssembler)
        
        [pb,assemblerTime] = updateDistribution(pb,distribution)
        
        pb = updateGreedySolver(pb,lhs,rhs)
        
        pb = updateModel(pb,model)
        
        pb = updateOperatorAssembler(pb,operatorAssembler)
        
        pb = updateUpdater(pb,lhs,rhs)
        
        [pb,time] = updateSource(pb,source)
        
    end
    
    methods (Access = private)
        
        pb = updateLocalSolver(pb,lhs,rhs)
        
    end
    
    methods (Static)
        
        pb = createRandom(varargin)
        
        function output = feOutputStructure(m,n)
            % output = feOutputStructure(m,n)
            if nargin < 2
                n = 1 ;
            end
            if numel(m) == 1
                sz = [m n] ;
            elseif numel(m) == 2
                sz = m ;
            else
                error('Bad input argument.')
            end
            fieldNames = {'error','residual','time','model','lhs','rhs',...
                'multiplier','operators'} ;
            fieldValues = repmat({[]},size(fieldNames)) ;
            fieldValues(1) = {repmat({[]},sz)} ;
            fields = [fieldNames ; fieldValues] ;
            output = struct(fields{:}) ;
            [output(:).operators] = deal(struct('diffOp',[],'srcOp',[],...
                'lhsBCOp',[],'rhsBCOp',[])) ;
        end
        
        function output = qpOutputStructure(sz,maxIter)
            % output = qpOutputStructure(sz,maxIter)
            if nargin < 2
                maxIter = 1 ;
            end
            if numel(sz)==1
                sz = [sz 1] ;
            end
            output = GreedyLinearSolver.outputStructure(sz,maxIter) ;
            [output.initTime] = deal([]) ;
        end
        
        function output = ddOutputStructure(maxIter)
            % output = ddOutputStructure(maxIter)
            fieldNames = {'flag','error','stagnation','time','initTime',...
                'iterTime','globalSolution','globalOutput',...
                'localSolution','localOutput'} ;
            arrayFields = [2 3 6] ;
            scalarFields = [1 4 5] ;
            cellFields = setdiff(1:numel(fieldNames),...
                union(arrayFields,scalarFields)) ;
            fieldValues = cell(size(fieldNames)) ;
            fieldValues(arrayFields) = {{zeros(maxIter,1)}} ;
            fieldValues(scalarFields) = {{[]}} ;
            fieldValues(cellFields) = {{cell(maxIter,1)}} ;
            fields = [fieldNames ; fieldValues] ;
            output = struct(fields{:}) ;
        end
        
    end
    
end