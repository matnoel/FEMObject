classdef QPProblem
    
    properties
        model
        K
        source
        bc
        tolerance
        penalty
        compress
        initialPoint
        useWeights
        maxIterationsGreedy
        maxIterationsUpdater
        maxIterationsLocal
        verboseGreedy
        verboseUpdater
        verboseLocal
    end
    
    methods( Access = public )
        %% Constructor methods
        function pb = QPProblem(model,conductivity,source,bc,varargin)
            % pb = QPProblem(model,conductivity,source,bc,varargin)
            if nargin < 3
                bc = QPBC(4) ; % Periodic BC
            end
            pb.model = model ;
            pb.K = conductivity ;
            pb.source = source ;
            if isnumeric(bc) && isscalar(bc)
                bc = QPBC(bc) ;
            end
            assert(isa(bc,'QPBC'),'Unknown BC type')
            pb.bc = bc ;
            p = ImprovedInputParser();
            addParameter(p,'tolerance',1e-3,@isnumeric);
            addParameter(p,'penalty',[],@isnumeric);
            addParameter(p,'compress',true,@islogical);
            addParameter(p,'initialPoint',[],@(x)isa(x,'TuckerLikeTensor'));
            addParameter(p,'useWeights',false,@islogical);
            addParameter(p,'maxIterationsGreedy',50,@isnumeric);
            addParameter(p,'maxIterationsUpdater',1,@isnumeric);
            addParameter(p,'maxIterationsLocal',10,@isnumeric);
            addParameter(p,'verboseGreedy',getVerbose(model),@islogical);
            addParameter(p,'verboseUpdater',false,@islogical);
            addParameter(p,'verboseLocal',false,@islogical);
            parse(p,varargin{:});
            pb = passMatchedArgsToProperties(p,pb);
        end
        
        function [sol,out,greedy] = solve(pb,lhs,rhs)
            if nargin < 3
                rhs = [] ;
                if nargin < 2
                    lhs = [] ;
                end
            end
            
            initClock = tic ;
            greedy = pb.greedySolver(lhs,rhs) ;
            x0 = pb.initialPoint ;
            if ~isempty(x0)
                initTime = toc(initClock) ;
                [x,out] = pb.mesoUpdate(greedy.A,greedy.b,x0) ;
                if ~out.flag
                    ifprint(pb.verboseUpdater,'Converged on initialisation\n');
                    sol = x;
                    out.initTime = initTime ;
                    return
                end
                if norm(x)>1e-16
                    greedy.x0 = x ;
                end
            end
            initTime = toc(initClock) ;
            
            [sol,out] = solve(greedy) ;
            
            out.initTime = initTime ;
        end
        
        function greedy = greedySolver(pb,lhs,rhs)
            if nargin < 3
                rhs = [] ;
                if nargin < 2
                    lhs = [] ;
                end
            end
            
            % Assemble what is necessary
            if isempty(lhs)
                if isempty(rhs)
                    [lhs,rhs] = swipOperator(pb.model,pb.K,pb.source,...
                        pb.bc,pb.penalty,pb.useWeights,pb.useWeights) ;
                else
                    lhs = swipOperator(pb.model,pb.K,pb.source,...
                        pb.bc,pb.penalty,pb.useWeights,pb.useWeights) ;
                end
            elseif isempty(rhs)
                [~,rhs] = swipOperator(pb.model,pb.K,pb.source,...
                    pb.bc,pb.penalty,pb.useWeights,pb.useWeights) ;
            end
            if pb.compress
                tr = Truncator('tolerance',getTolSVD(pb.model)) ;
                lhs = tr.truncate(lhs) ;
                rhs = tr.truncate(rhs) ;
            end
            
            % Create greedy solver
            tol = pb.tolerance ;
            maxIterG = min(getCellNb(pb.model),pb.maxIterationsGreedy) ;
            up = pb.updater(lhs,rhs) ;
            solver = pb.localSolver(lhs,rhs) ;
            greedy = GreedyLinearSolver(lhs,rhs,...
                'localSolver',solver,...
                'update',@(x)updateTSpaceByALS(up,x),...
                'x0',[],...
                'initialUpdate',false,...
                'tolerance',tol, ...
                'maxIterations',maxIterG,...
                'stagnation',min(1e-6,tol/100),...
                'minimizeResidual',false,...
                'checkResidual',1,... % Set to N to check residual every N steps
                'display',pb.verboseGreedy);
        end
        
        function solver = localSolver(pb,lhs,rhs)
            solver = RankOneALSLinearSolver(lhs,rhs, ...
                'maxIterations',10,...
                'stagnation',min(1e-6,pb.tolerance/100),...
                'display',pb.verboseLocal,...
                'x0',TuckerLikeTensor.ones(pb.model.tensorSize),...
                'tolerance',1e-3) ; % tolerance is unused
        end
        
        function updater = updater(pb,lhs,rhs)
            tol = pb.tolerance ;
            updater = TuckerLikeTensorALSLinearSolver(lhs,rhs,...
                'tolerance',tol,...
                'maxIterationsTSpace',pb.maxIterationsUpdater,...
                'stagnationTSpace',min(1e-6,tol/100),...
                'display',pb.verboseUpdater) ;
        end
        
        function [x,out] = mesoUpdate(pb,lhs,rhs,x0)
            % [sol,out] = mesoUpdate(pb,lhs,rhs,initialPoint)
            if nargin < 4 || isempty(x0)
                x0 = pb.initialPoint ;
            end
            
            assert(~isempty(x0),'Initial point is empty')
            
            clock0 = tic;
            out = GreedyLinearSolver.outputStructure(1,1) ;
            updater = pb.updater(lhs,rhs) ;
            
            xAx = dotWithMetrics(x0.space,x0.space,updater.A.space);
            xb  = dot(updater.b.space,x0.space);
            x = updateTSpaceDim(updater,x0,1,xAx,xb);
            if x.order == 3
                xAx(1) = dotWithMetrics(x.space,x.space,s.A.space,1);
                xb(1) = dot(s.b.space,x.space,1);
                x = updateTSpaceDim(updater,x,2,xAx,xb);
            end
            
            out.error = norm(updater.A*x-updater.b)/norm(updater.b);
            out.errvec = out.error ;
            % no stagnation criterion for initial update
            out.flag = out.error >= updater.tolerance || isnan(out.error);
            ifprint(pb.verboseUpdater,'mesoUpdate: error %.2d\n',out.error)
            out.iter = 0 ;
            out.time = toc(clock0);
            out.times = out.time ;
        end
        
        %% External methods signatures
        
        [feModel,mesherTime] = buildFEModel(pb,method)
        
        K = getCompleteK(pb,fullModel)
        
        source = getCompleteSource(pb,fullModel)
                
        [Khomo,output] = homogenise(pb,K,stagnationTol,standardDeviationTol,method)
        
        [solution,output] = solveFE(pb,useMultiplier,meshMethod)
        
    end
    
    methods (Static)
       
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
        
    end
    
end