classdef QPDiffusionProblemTest < matlab.unittest.TestCase
    
    properties
        testTol
        testCrit
    end
    
    methods (Test)
        
        function construction(testCase)
            noError = true ;
            for iter = 1:10
                try
                    QPDiffusionProblem.createRandom('verbose',false) ;
                catch
                    noError = false ;
                    break
                end
            end
            testCase.verifyTrue(noError) ;
        end
        
        %% QP2 vs FEM
        
        % % BC
        % Dirichlet
        
        function QP2vsFEMDirichletHomogeneous(testCase)
            model = QPDiffusionProblemTest.testModel(2,[1 1]) ;
            cA = QPDiffusionProblemTest.testCA(model,1,[]) ;
            bc = QPDiffusionProblemTest.testBC(model,1,2.5,1.1e3) ;
            dA = QPDiffusionProblemTest.testDA(cA,{bc},{},0) ;
            pb = QPDiffusionProblem('operatorAssembler',dA,...
                'tolerance',testCase.testTol) ;
            compareQPvsFEM(testCase,pb) ;
        end
        function QP2vsFEMDirichletHeterogeneous(testCase)
            model = QPDiffusionProblemTest.testModel(2,[3 2]) ;
            cA = QPDiffusionProblemTest.testCA(model,2,5) ;
            bc = QPDiffusionProblemTest.testBC(model,1,2.5,1.1e3) ;
            dA = QPDiffusionProblemTest.testDA(cA,{bc},{},0) ;
            pb = QPDiffusionProblem('operatorAssembler',dA,...
                'tolerance',testCase.testTol) ;
            compareQPvsFEM(testCase,pb) ;
        end
        
        % Neumann
        
        function QP2vsFEMNeumannHomogeneous(testCase)
            model = QPDiffusionProblemTest.testModel(2,[1 1]) ;
            cA = QPDiffusionProblemTest.testCA(model,1,[]) ;
            bc = QPDiffusionProblemTest.testBC(model,2,2.5,[]) ;
            dA = QPDiffusionProblemTest.testDA(cA,{bc},{},0) ;
            pb = QPDiffusionProblem('operatorAssembler',dA,...
                'tolerance',testCase.testTol) ;
            compareQPvsFEM(testCase,pb) ;
        end
        function QP2vsFEMNeumannHeterogeneous(testCase)
            model = QPDiffusionProblemTest.testModel(2,[3 2]) ;
            cA = QPDiffusionProblemTest.testCA(model,2,5) ;
            bc = QPDiffusionProblemTest.testBC(model,2,2.5,[]) ;
            dA = QPDiffusionProblemTest.testDA(cA,{bc},{},0) ;
            pb = QPDiffusionProblem('operatorAssembler',dA,...
                'tolerance',testCase.testTol) ;
            compareQPvsFEM(testCase,pb) ;
        end
        
        % Robin
        
        function QP2vsFEMRobinHomogeneous(testCase)
            model = QPDiffusionProblemTest.testModel(2,[1 1]) ;
            cA = QPDiffusionProblemTest.testCA(model,1,[]) ;
            bc = QPDiffusionProblemTest.testBC(model,3,2.5,1.1) ;
            dA = QPDiffusionProblemTest.testDA(cA,{bc},{},0) ;
            pb = QPDiffusionProblem('operatorAssembler',dA,...
                'tolerance',testCase.testTol) ;
            compareQPvsFEM(testCase,pb) ;
        end
        function QP2vsFEMRobinHeterogeneous(testCase)
            model = QPDiffusionProblemTest.testModel(2,[3 2]) ;
            cA = QPDiffusionProblemTest.testCA(model,2,5) ;
            bc = QPDiffusionProblemTest.testBC(model,3,2.5,1.1) ;
            dA = QPDiffusionProblemTest.testDA(cA,{bc},{},0) ;
            pb = QPDiffusionProblem('operatorAssembler',dA,...
                'tolerance',testCase.testTol) ;
            compareQPvsFEM(testCase,pb) ;
        end
        
        % "Cauchy"
        
%         function QP2vsFEMCauchyHomogeneous(testCase)
%             model = QPDiffusionProblemTest.testModel(2,[1 1]) ;
%             cA = QPDiffusionProblemTest.testCA(model,1,[]) ;
%             bc = QPDiffusionProblemTest.testBC(model,4,{2.5 5.5},1.1e3) ;
%             dA = QPDiffusionProblemTest.testDA(cA,{bc},{},0) ;
%             pb = QPDiffusionProblem('operatorAssembler',dA,...
%                 'tolerance',testCase.testTol) ;
%             compareQPvsFEM(testCase,pb) ;
%         end
%         function QP2vsFEMCauchyHeterogeneous(testCase)
%             model = QPDiffusionProblemTest.testModel(2,[3 2]) ;
%             cA = QPDiffusionProblemTest.testCA(model,2,5) ;
%             bc = QPDiffusionProblemTest.testBC(model,4,{2.5 5.5},1.1e3) ;
%             dA = QPDiffusionProblemTest.testDA(cA,{bc},{},0) ;
%             pb = QPDiffusionProblem('operatorAssembler',dA,...
%                 'tolerance',testCase.testTol) ;
%             compareQPvsFEM(testCase,pb) ;
%         end  
        
        % Periodic
        
        function QP2vsFEMPeriodicHomogeneous(testCase)
            model = QPDiffusionProblemTest.testModel(2,[1 1]) ;
            cA = QPDiffusionProblemTest.testCA(model,1,[]) ;
            bc = QPDiffusionProblemTest.testBC(model,5,[],[]) ;
            dA = QPDiffusionProblemTest.testDA(cA,{bc},{},6.7) ;
            pb = QPDiffusionProblem('operatorAssembler',dA,...
                'tolerance',testCase.testTol) ;
            compareQPvsFEM(testCase,pb) ;
        end
        function QP2vsFEMPeriodicHeterogeneous(testCase)
            model = QPDiffusionProblemTest.testModel(2,[3 2]) ;
            cA = QPDiffusionProblemTest.testCA(model,2,5) ;
            bc = QPDiffusionProblemTest.testBC(model,5,[],[]) ;
            dA = QPDiffusionProblemTest.testDA(cA,{bc},{},6.7) ;
            pb = QPDiffusionProblem('operatorAssembler',dA,...
                'tolerance',testCase.testTol) ;
            compareQPvsFEM(testCase,pb) ;
        end  
        
        % % Source
        
        function QP2vsFEMPeriodicExpSrc(testCase)
            model = QPDiffusionProblemTest.testModel(2,[3 2]) ;
            cA = QPDiffusionProblemTest.testCA(model,2,5) ;
            bc = QPDiffusionProblemTest.testBC(model,5,[],[]) ;
            src = QPDiffusionProblemTest.testSource(model) ;
            dA = QPDiffusionProblemTest.testDA(cA,{bc},{},src) ;
            pb = QPDiffusionProblem('operatorAssembler',dA,...
                'tolerance',testCase.testTol) ;
            compareQPvsFEM(testCase,pb) ;
        end
        
        % % Patches
        % SWIP
                
%         function QP2vsFEMSWIPPatchQPSolver(testCase) % with PBC
%             model = QPDiffusionProblemTest.testModel(2,[4 3]) ;
%             cA = QPDiffusionProblemTest.testCA(model,1,[]) ;
%             bc = QPDiffusionProblemTest.testBC(model,2,2.5,[]) ;
%             cAp = QPDiffusionProblemTest.testCA(model,2,[6;7]) ;
%             ps = QPDiffusionProblemTest.testPatches(cAp,1,{{[6;7]}},1) ;
%             dA = QPDiffusionProblemTest.testDA(cA,{bc},ps,1.8) ;
%             pb = QPDiffusionProblem('operatorAssembler',dA,...
%                 'tolerance',testCase.testTol) ;
%             compareQPvsFEM(testCase,pb) ;
%         end
                
%         function QP2vsFEMSWIPPatchFEMSolver(testCase)
%             model = QPDiffusionProblemTest.testModel(2,[4 3]) ;
%             cA = QPDiffusionProblemTest.testCA(model,1,[]) ;
%             bc = QPDiffusionProblemTest.testBC(model,2,2.5,[]) ;
%             cAp = QPDiffusionProblemTest.testCA(model,2,[6;7]) ;
%             ps = QPDiffusionProblemTest.testPatches(cAp,1,{{[6;7]}},2) ;
%             dA = QPDiffusionProblemTest.testDA(cA,{bc},ps,1.8) ;
%             pb = QPDiffusionProblem('operatorAssembler',dA,...
%                 'tolerance',testCase.testTol) ;
%             compareQPvsFEM(testCase,pb) ;
%         end
        
        % Neumann-Dirichlet
                
%         function QP2vsFEMNDPatchQPSolver(testCase)
%             model = QPDiffusionProblemTest.testModel(2,[4 3]) ;
%             cA = QPDiffusionProblemTest.testCA(model,1,[]) ;
%             bc = QPDiffusionProblemTest.testBC(model,1,0,[]) ;
%             cAp = QPDiffusionProblemTest.testCA(model,2,[6;7]) ;
%             ps = QPDiffusionProblemTest.testPatches(cAp,2,{{[6;7]}},1) ;
%             dA = QPDiffusionProblemTest.testDA(cA,{bc},ps,1.8) ;
%             pb = QPDiffusionProblem('operatorAssembler',dA,...
%                 'tolerance',testCase.testTol) ;
%             compareQPvsFEM(testCase,pb) ;
%         end
                
%         function QP2vsFEMNDPatchFEMSolver(testCase)
%             model = QPDiffusionProblemTest.testModel(2,[4 3]) ;
%             cA = QPDiffusionProblemTest.testCA(model,1,[]) ;
%             bc = QPDiffusionProblemTest.testBC(model,1,0,[]) ;
%             cAp = QPDiffusionProblemTest.testCA(model,2,[6;7]) ;
%             ps = QPDiffusionProblemTest.testPatches(cAp,2,{{[6;7]}},2) ;
%             dA = QPDiffusionProblemTest.testDA(cA,{bc},ps,1.8) ;
%             pb = QPDiffusionProblem('operatorAssembler',dA,...
%                 'tolerance',testCase.testTol) ;
%             compareQPvsFEM(testCase,pb) ;
%         end
        
        % Lagrange
                
%         function QP2vsFEMLagrangePatchFEMSolver(testCase)
%             model = QPDiffusionProblemTest.testModel(2,[4 3]) ;
%             cA = QPDiffusionProblemTest.testCA(model,1,[]) ;
%             bc = QPDiffusionProblemTest.testBC(model,2,2.5,[]) ;
%             cAp = QPDiffusionProblemTest.testCA(model,2,[6;7]) ;
%             ps = QPDiffusionProblemTest.testPatches(cAp,3,{{[6;7]}},2) ;
%             dA = QPDiffusionProblemTest.testDA(cA,{bc},ps,1.8) ;
%             pb = QPDiffusionProblem('operatorAssembler',dA,...
%                 'tolerance',testCase.testTol) ;
%             compareQPvsFEM(testCase,pb) ;
%         end
        
        
        
        %% QP3 vs FEM
        
        % % BC
        % Dirichlet
        
        function QP3vsFEMDirichletHomogeneous(testCase)
            model = QPDiffusionProblemTest.testModel(3,[1 1]) ;
            cA = QPDiffusionProblemTest.testCA(model,1,[]) ;
            bc = QPDiffusionProblemTest.testBC(model,1,2.5,1.1e3) ;
            dA = QPDiffusionProblemTest.testDA(cA,{bc},{},0) ;
            pb = QPDiffusionProblem('operatorAssembler',dA,...
                'tolerance',testCase.testTol) ;
            compareQPvsFEM(testCase,pb) ;
        end
        function QP3vsFEMDirichletHeterogeneous(testCase)
            model = QPDiffusionProblemTest.testModel(3,[3 2]) ;
            cA = QPDiffusionProblemTest.testCA(model,2,5) ;
            bc = QPDiffusionProblemTest.testBC(model,1,2.5,1.1e3) ;
            dA = QPDiffusionProblemTest.testDA(cA,{bc},{},0) ;
            pb = QPDiffusionProblem('operatorAssembler',dA,...
                'tolerance',testCase.testTol) ;
            compareQPvsFEM(testCase,pb) ;
        end
        
        % Neumann
        
        function QP3vsFEMNeumannHomogeneous(testCase)
            model = QPDiffusionProblemTest.testModel(3,[1 1]) ;
            cA = QPDiffusionProblemTest.testCA(model,1,[]) ;
            bc = QPDiffusionProblemTest.testBC(model,2,2.5,[]) ;
            dA = QPDiffusionProblemTest.testDA(cA,{bc},{},0) ;
            pb = QPDiffusionProblem('operatorAssembler',dA,...
                'tolerance',testCase.testTol) ;
            compareQPvsFEM(testCase,pb) ;
        end
        function QP3vsFEMNeumannHeterogeneous(testCase)
            model = QPDiffusionProblemTest.testModel(3,[3 2]) ;
            cA = QPDiffusionProblemTest.testCA(model,2,5) ;
            bc = QPDiffusionProblemTest.testBC(model,2,2.5,[]) ;
            dA = QPDiffusionProblemTest.testDA(cA,{bc},{},0) ;
            pb = QPDiffusionProblem('operatorAssembler',dA,...
                'tolerance',testCase.testTol) ;
            compareQPvsFEM(testCase,pb) ;
        end
        
        % Robin
        
        function QP3vsFEMRobinHomogeneous(testCase)
            model = QPDiffusionProblemTest.testModel(3,[1 1]) ;
            cA = QPDiffusionProblemTest.testCA(model,1,[]) ;
            bc = QPDiffusionProblemTest.testBC(model,3,2.5,1.1) ;
            dA = QPDiffusionProblemTest.testDA(cA,{bc},{},0) ;
            pb = QPDiffusionProblem('operatorAssembler',dA,...
                'tolerance',testCase.testTol) ;
            compareQPvsFEM(testCase,pb) ;
        end
        function QP3vsFEMRobinHeterogeneous(testCase)
            model = QPDiffusionProblemTest.testModel(3,[3 2]) ;
            cA = QPDiffusionProblemTest.testCA(model,2,5) ;
            bc = QPDiffusionProblemTest.testBC(model,3,2.5,1.1) ;
            dA = QPDiffusionProblemTest.testDA(cA,{bc},{},0) ;
            pb = QPDiffusionProblem('operatorAssembler',dA,...
                'tolerance',testCase.testTol) ;
            compareQPvsFEM(testCase,pb) ;
        end
        
        % "Cauchy"
        
%         function QP3vsFEMCauchyHomogeneous(testCase)
%             model = QPDiffusionProblemTest.testModel(3,[1 1]) ;
%             cA = QPDiffusionProblemTest.testCA(model,1,[]) ;
%             bc = QPDiffusionProblemTest.testBC(model,4,{2.5 5.5},1.1e3) ;
%             dA = QPDiffusionProblemTest.testDA(cA,{bc},{},0) ;
%             pb = QPDiffusionProblem('operatorAssembler',dA,...
%                 'tolerance',testCase.testTol) ;
%             compareQPvsFEM(testCase,pb) ;
%         end
%         function QP3vsFEMCauchyHeterogeneous(testCase)
%             model = QPDiffusionProblemTest.testModel(3,[3 2]) ;
%             cA = QPDiffusionProblemTest.testCA(model,2,5) ;
%             bc = QPDiffusionProblemTest.testBC(model,4,{2.5 5.5},1.1e3) ;
%             dA = QPDiffusionProblemTest.testDA(cA,{bc},{},0) ;
%             pb = QPDiffusionProblem('operatorAssembler',dA,...
%                 'tolerance',testCase.testTol) ;
%             compareQPvsFEM(testCase,pb) ;
%         end  
        
        % Periodic
        
        function QP3vsFEMPeriodicHomogeneous(testCase)
            model = QPDiffusionProblemTest.testModel(3,[1 1]) ;
            cA = QPDiffusionProblemTest.testCA(model,1,[]) ;
            bc = QPDiffusionProblemTest.testBC(model,5,[],[]) ;
            dA = QPDiffusionProblemTest.testDA(cA,{bc},{},6.7) ;
            pb = QPDiffusionProblem('operatorAssembler',dA,...
                'tolerance',testCase.testTol) ;
            compareQPvsFEM(testCase,pb) ;
        end
        function QP3vsFEMPeriodicHeterogeneous(testCase)
            model = QPDiffusionProblemTest.testModel(3,[3 2]) ;
            cA = QPDiffusionProblemTest.testCA(model,2,5) ;
            bc = QPDiffusionProblemTest.testBC(model,5,[],[]) ;
            dA = QPDiffusionProblemTest.testDA(cA,{bc},{},6.7) ;
            pb = QPDiffusionProblem('operatorAssembler',dA,...
                'tolerance',testCase.testTol) ;
            compareQPvsFEM(testCase,pb) ;
        end  
        
        % % Source
        
        function QP3vsFEMPeriodicExpSrc(testCase)
            model = QPDiffusionProblemTest.testModel(3,[3 2]) ;
            cA = QPDiffusionProblemTest.testCA(model,2,5) ;
            bc = QPDiffusionProblemTest.testBC(model,5,[],[]) ;
            src = QPDiffusionProblemTest.testSource(model) ;
            dA = QPDiffusionProblemTest.testDA(cA,{bc},{},src) ;
            pb = QPDiffusionProblem('operatorAssembler',dA,...
                'tolerance',testCase.testTol) ;
            compareQPvsFEM(testCase,pb) ;
        end
        
        % % Patches
        % SWIP
                
%         function QP3vsFEMSWIPPatchQPSolver(testCase) % with PBC
%             model = QPDiffusionProblemTest.testModel(3,[4 3]) ;
%             cA = QPDiffusionProblemTest.testCA(model,1,[]) ;
%             bc = QPDiffusionProblemTest.testBC(model,2,2.5,[]) ;
%             cAp = QPDiffusionProblemTest.testCA(model,2,[6;7]) ;
%             ps = QPDiffusionProblemTest.testPatches(cAp,1,{{[6;7]}},1) ;
%             dA = QPDiffusionProblemTest.testDA(cA,{bc},ps,1.8) ;
%             pb = QPDiffusionProblem('operatorAssembler',dA,...
%                 'tolerance',testCase.testTol) ;
%             compareQPvsFEM(testCase,pb) ;
%         end
                 
%         function QP3vsFEMSWIPPatchFEMSolver(testCase)
%             model = QPDiffusionProblemTest.testModel(3,[4 3]) ;
%             cA = QPDiffusionProblemTest.testCA(model,1,[]) ;
%             bc = QPDiffusionProblemTest.testBC(model,2,2.5,[]) ;
%             cAp = QPDiffusionProblemTest.testCA(model,2,[6;7]) ;
%             ps = QPDiffusionProblemTest.testPatches(cAp,1,{{[6;7]}},2) ;
%             dA = QPDiffusionProblemTest.testDA(cA,{bc},ps,1.8) ;
%             pb = QPDiffusionProblem('operatorAssembler',dA,...
%                 'tolerance',testCase.testTol) ;
%             compareQPvsFEM(testCase,pb) ;
%         end
        
        % Neumann-Dirichlet
                
%         function QP3vsFEMNDPatchQPSolver(testCase)
%             model = QPDiffusionProblemTest.testModel(3,[4 3]) ;
%             cA = QPDiffusionProblemTest.testCA(model,1,[]) ;
%             bc = QPDiffusionProblemTest.testBC(model,2,2.5,[]) ;
%             cAp = QPDiffusionProblemTest.testCA(model,2,[6;7]) ;
%             ps = QPDiffusionProblemTest.testPatches(cAp,2,{{[6;7]}},1) ;
%             dA = QPDiffusionProblemTest.testDA(cA,{bc},ps,1.8) ;
%             pb = QPDiffusionProblem('operatorAssembler',dA,...
%                 'tolerance',testCase.testTol) ;
%             compareQPvsFEM(testCase,pb) ;
%         end
                
%         function QP3vsFEMNDPatchFEMSolver(testCase)
%             model = QPDiffusionProblemTest.testModel(3,[4 3]) ;
%             cA = QPDiffusionProblemTest.testCA(model,1,[]) ;
%             bc = QPDiffusionProblemTest.testBC(model,2,2.5,[]) ;
%             cAp = QPDiffusionProblemTest.testCA(model,2,[6;7]) ;
%             ps = QPDiffusionProblemTest.testPatches(cAp,2,{{[6;7]}},2) ;
%             dA = QPDiffusionProblemTest.testDA(cA,{bc},ps,1.8) ;
%             pb = QPDiffusionProblem('operatorAssembler',dA,...
%                 'tolerance',testCase.testTol) ;
%             compareQPvsFEM(testCase,pb) ;
%         end
        
        % Lagrange
                
%         function QP3vsFEMLagrangePatchFEMSolver(testCase)
%             model = QPDiffusionProblemTest.testModel(3,[4 3]) ;
%             cA = QPDiffusionProblemTest.testCA(model,1,[]) ;
%             bc = QPDiffusionProblemTest.testBC(model,2,2.5,[]) ;
%             cAp = QPDiffusionProblemTest.testCA(model,2,[6;7]) ;
%             ps = QPDiffusionProblemTest.testPatches(cAp,3,{{[6;7]}},2) ;
%             dA = QPDiffusionProblemTest.testDA(cA,{bc},ps,1.8) ;
%             pb = QPDiffusionProblem('operatorAssembler',dA,...
%                 'tolerance',testCase.testTol) ;
%             compareQPvsFEM(testCase,pb) ;
%         end
         
    end
    
    %% Static methods
    methods (Static)
        
        function model = testModel(order,cellNum)
            model = QPModel('order',order,'cellNum',cellNum,'verbose',...
                false,'cellSize',[1 1],'elementSize',[1 1]/10,...
                'tolSVD',1e-6) ;
        end
        
        function cA = testCA(model,type,cells)
            switch type
                case 1 % Homogeneous
                    cA = QPConductivityAssembler.homogeneous(1.5,model) ;
                case 2 % Discs
                    pattern = struct('name',{'uniform','disc'},'value',...
                        {1 8},'size',{[],.2},'center',{[],[.5 .5]},...
                        'offset',{[],[]}) ;
                    cellNb = getCellNb(model) ;
                    dist = {setdiff(1:cellNb,cells)' cells} ;
                    cA = QPConductivityAssembler('model',model,...
                        'patterns',pattern,'patternsTable',[1 1 ; 0 1],...
                        'distribution',dist) ;
            end
            cA = assemble(cA) ;
        end
        
        function src = testSource(model)
            dsz = getDomainSize(model) ;
            src = @(x) exp(-10*((x(:,1)-dsz(1)/2).^2+(x(:,2)-dsz(2)/2).^2)) ;
        end
        
        function bc = testBC(model,type,value,penalty)
            bc = QPBC('model',model,'type',type,'value',value,'penalty',...
                penalty,'cells',boundaryCells(model),'isSIP',false,...
                'useAverageWeights',false,'useStabilisationWeights',false) ;
        end
        
        function dA = testDA(cA,bc,ps,source)
            dA = QPDiffusionAssembler('conductivityAssembler',cA,'bc',bc,...
                'source',source,'penalty',[],'patches',ps,'useCompression',...
                false,'constantNullification','full','useAverageWeights',...
                true,'useStabilisationWeights',true);
            dA = assemble(dA) ;
        end
        
        function ps = testPatches(cA,type,cells,solver)
            cells = formatIndex(getOrder(cA),getCellNum(cA),cells) ;
            Kpatch = QPDiffusionProblemTest.subConductivity(getConductivity(cA),cells) ;
            ps = QPPatches('conductivityAssembler',cA,'cells',cells,...
                'useCompression',false,'type',type,'solver',solver,...
                'conductivity',Kpatch) ;
        end
        
        function Kpatch = subConductivity(K,pCells)
            Kpatch = cell(numel(pCells),1) ;
            for p = 1:numel(pCells)
                pSz = size(pCells{p}{1}) ;
                subs = [mat2cell(pCells{p}{1},pSz(1),ones(1,pSz(2))), {':'}] ;
                subs = cellfun(@unique,subs,'UniformOutput',false) ;
                Kpatch{p} = subTensor(K,subs{:}) ;
            end
        end
        
    end
    
    %%    
    methods
        
        function pbT = QPDiffusionProblemTest(tol,errFun)
            if nargin < 2
%                 errFun = 'residual' ;
                errFun = @(a,b) max(-1,norm(full(b-a))/norm(full(a))) ;
                if nargin < 1
                    tol = 1e-4 ;                    
                end
            end
            pbT.testTol = tol ;
            pbT.testCrit = errFun ;
        end
        
        function [] = compareQPvsFEM(testCase,pb)
            [solQP,~,pb] = solve(pb) ;
            [solFE,outFE] = solveFEM(pb) ;
            errorQP = compare2FE(pb,solQP,solFE,testCase.testCrit,...
                'feModel',outFE.model) ;
            %/DEBUG
%             disp(errorQP)
%             if isempty(getPatches(pb))
%                 model = getModel(pb) ;
%                 feCoord = getcoord(getnode(outFE.model)) ;
%                 qpCoord = smooth(model,getDomainCoord(model)) ;
%                 solFE = relativeSort(solFE,feCoord,qpCoord) ;                
%                 figure
%                 subplot(1,2,1)
%                 plot(model,solQP)
%                 subplot(1,2,2)
%                 plot(model,solFE)
%             end
            %DEBUG\
            testCase.verifyLessThan(errorQP,testCase.testTol) ;
        end
        
    end
end