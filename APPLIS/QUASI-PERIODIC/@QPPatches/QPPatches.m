classdef QPPatches
    
    properties( SetAccess = public, GetAccess = public )
        patches
        cells
        conductivityAssembler
        useCompression
        lhsOperator
        rhsOperator
    end
    
    methods( Access = public )
   
        % Constructor
        function ps = QPPatches(varargin)
            ip = ImprovedInputParser();
            addParameter(ip,'patches',{});
            addParameter(ip,'conductivityAssembler',[]);
            addParameter(ip,'cells',{});
            addParameter(ip,'useCompression',false);
            addParameter(ip,'lhsOperator',{});
            addParameter(ip,'rhsOperator',{});            
            if ~verLessThan('matlab','8.2') % compatibility (>R2013a)
                ip.PartialMatching = false ;
            end
            parse(ip,varargin{:});
            ps = passMatchedArgsToProperties(ip,ps);
            ps = setCells(ps) ;
            pCells = getCells(ps) ;
            assert(checkCells(ps)>0,'Patches cells incorrect.')
            if ischarin('type',varargin) % Create new patches
                pCells(1:numel(ps)) = [] ; % remove cells related to patches
                varargin = setcharin('cells',varargin,{pCells}) ;
                varargin = setcharin('model',varargin,getModel(ps)) ;
                K = getcharin('conductivity',varargin) ;
                if ~iscell(K) % Ensure cell array of TuckerLikeTensors
                    K = repmat({K},numel(pCells),1) ;
                    varargin = setcharin('conductivity',varargin,{K}) ;
                end
                newPatches = QPPatch.batch(varargin{:}) ;
                ps = setPatches(ps,[getPatches(ps) ; newPatches]) ;
            end
            realCA = realConductivityAssembler(ps,getConductivityAssembler(ps),1) ;
            ps = setConductivityAssembler(ps,realCA) ;
        end
        
        %% Get methods
        
        % Primary
        
        function patchesCellArray = getPatches(Ps)
            patchesCellArray = Ps.patches ;
        end
        
        function cells = getCells(ps,patchNum)
            if nargin < 2
                cells = ps.cells ;
            else
                cells = ps.cells ;
                cells = cells(patchNum) ;
            end
        end
        
        function cA = getConductivityAssembler(ps)
            cA = ps.conductivityAssembler ;
        end
        
        function lhsOp = getLHSOperator(ps)
            lhsOp = ps.lhsOperator ;
        end
        
        function rhsOp = getRHSOperator(ps)
            rhsOp = ps.rhsOperator ;
        end
        
        function useCompression = getUseCompression(ps)
            useCompression = ps.useCompression ;
        end
        
        % Secondary
        
        function patchesNb = numel(Ps)
            patchesNb = numel(getPatches(Ps)) ;
        end
        
        function flag = isempty(ps,patchNum)
           if nargin == 1
               patchNum = 1:numel(ps) ;
           end
            flag = true ;
            pCells = getCells(ps) ;
            for i = patchNum
                if iscell(pCells{i})
                    flag = flag && all(cellfun(@isempty,pCells{i})) ;
                else
                    flag = flag && isempty(pCells{i}) ;
                end
            end
        end
        
        function patch = getPatch(Ps,patchNum)
            patch = getPatches(Ps) ;
            if numel(patchNum)==1
                patch = patch{patchNum} ;
            else
                patch = patch(patchNum) ;
            end
        end
        
        function cellList = getCellList(ps)
           pcells = getCells(ps) ;
           cellList = [] ;
           for c = 1:numel(pcells)
               cellList = [cellList ; cell2mat(pcells{c}(:))] ;
           end
        end
        
        % References to QPPatch
        
        function type = getType(ps,patchNum)
            if nargin < 2
                patchNum = 1:numel(ps) ;
            end
            type = cell(numel(patchNum),1) ;
            p = getPatch(ps,patchNum) ;
            for i = patchNum
                type{i} = getType(p{i}) ;
            end
            if numel(type)==1
                type = type{1} ;
            end
        end
        
        function pModel = getPatchModel(ps,patchNum)
            if nargin < 2
                patchNum = 1:numel(ps) ;
            end
            pModel = cell(numel(patchNum),1) ;
            p = getPatch(ps,patchNum) ;
            if numel(patchNum)==1
                pModel = getModel(p) ;
            else
                for i = patchNum
                    pModel{i} = getModel(p{i}) ;
                end
            end
        end
        
        function solver = getSolver(ps,patchNum)
            if nargin < 2
                patchNum = 1:numel(ps) ;
            end
            solver = cellfun(@getSolver,getPatch(ps,patchNum)) ;
        end
        
        function conductivity = getPatchConductivity(ps,patchNum)
            if nargin < 2
                patchNum = 1:numel(ps) ;
            end
            conductivity = cellfun(@getConductivity,getPatch(ps,patchNum),...
                'UniformOutput',false) ;
            if numel(conductivity)==1
                conductivity = conductivity{1} ;
            end
        end
        
        % References to QPConductivityAssembler
        
        function model = getModel(ps)
            model = getModel(getConductivityAssembler(ps)) ;
        end
        
        function bounds = getConductivityBounds(ps,choice)
            if nargin==1 ; choice = 'all' ; end
            bounds = getConductivityBounds(getConductivityAssembler(ps),choice) ;
        end
        
        function conductivity = getConductivity(ps)
            conductivity = getConductivity(getConductivityAssembler(ps)) ;
        end
        
        % References to QPModel
        
        function order = getOrder(ps)
            order = getOrder(getModel(ps)) ;
        end
        
        function cellNum = getCellNum(ps)
            cellNum = getCellNum(getModel(ps)) ;
        end
        
        function verbose = getVerbose(ps)
            verbose = getVerbose(getModel(ps)) ;
        end
        
        %% Set methods
        
        % Primary
        
        function ps = setPatches(ps,newPatches)
            ps.patches = newPatches ;
        end
        
        function ps = setCells(ps,cells)
            if nargin < 2
                cells = formatCells(ps,getCells(ps)) ;
            end
            ps.cells = cells ;
        end
        
        function ps = setConductivityAssembler(ps,cA)
            ps.conductivityAssembler = cA ;
        end
        
        function ps = setLHSOperator(ps,lhsOp)
            ps.lhsOperator = lhsOp ;
        end
        
        function ps = setRHSOperator(ps,rhsOp)
            ps.rhsOperator = rhsOp ;
        end
        
        function ps = setUseCompression(ps,newUseCompression)
            ps.useCompression = newUseCompression ;
        end
        
        % Secondary
        
        function ps = addPatch(ps,patch,pcells)
            if ~iscell(pcells{1})
                pcells = {pcells} ;
            end
            if ~iscell(patch)
                patch = {patch} ;
            end
            ps = setCells(ps,[getCells(ps);pcells]) ;
            ps = setPatches(ps,[getPatches(ps);patch]) ;
        end
        
        function ps = removePatch(ps,patchNum)
            if nargin == 1
                patchNum = 1:numel(ps) ;
            end
            patchesLeft = setdiff(1:numel(ps),patchNum) ;
            ps = setPatches(ps,getPatch(ps,patchesLeft)) ;
            ps = setCells(ps,getCells(ps,patchesLeft)) ;
        end
        
        function ps = addCells(ps,newCells,patchNum)
            assert(numel(patchNum)==1 || ...
                (iscell(newCells)&&numel(newCells)==numel(patchNum)),...
                'Incorrect input argument') ;
            if ~iscell(newCells)
                newCells = {newCells} ;
            end
            pCells = getCells(ps) ;
            for i = 1:numel(patchNum)
                pCells(patchNum(i)) = [pCells(patchNum(i)) ; newCells(i)] ;
            end
            ps = setCells(ps,pCells) ;
        end
        
        function ps = setPatchConductivity(ps,conductivity,patchNum)
            if nargin < 3
                patchNum = 1:numel(ps) ;
            end
            if ~iscell(conductivity)
                conductivity = repmat({conductivity},numel(patchNum),1) ;
            end
            patches = getPatches(ps) ;
            for i = 1:numel(patchNum)
                patches{patchNum(i)} = setConductivity(patches{patchNum(i)},...
                    conductivity{i}) ;
            end
            ps = setPatches(ps,patches) ;
        end
        
        % References
        
        function ps = setVerbose(ps,newVerbose)
            cA = setVerbose(getConductivityAssembler(ps),newVerbose) ;
            ps = setConductivityAssembler(ps,cA) ;
        end
        
        function ps = setModel(ps,newModel)
            cA = setModel(getConductivityAssembler(ps),newModel) ;
            ps = setConductivityAssembler(ps,cA) ;
        end
        
        %% External methods signatures
        
        [ps,time] = assemble(ps,microOperators,dGPenalty)
        
        flag = checkCells(ps)
        
        flag = checkOverlap(ps)
        
        cells = formatCells(ps,cells) ;
        
        realCA = realConductivityAssembler(ps,globalCA,compress)
        
        [sol,value,flux,output] = solve(ps,microOp,globalSource,globalSolution,tolerance,penalty)
        
        ps = updateConductivityAssembler(ps,globalCA)
        
        ps = updateModel(ps,model)
       
    end
    
end