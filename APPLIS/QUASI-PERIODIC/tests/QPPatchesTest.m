classdef QPPatchesTest < matlab.unittest.TestCase
    
    properties
        testTol = 1e-9 ;
        testCrit = @(a,b) max(-1,normest(b-a)/normest(a)) ;
    end
    
    methods (Test)
        
%         function QP3vsQP2SWIP1x1Homogeneous(testCase)
%             model2 = QPPatchesTest.testModel(2,[3 3]) ;
%             cA2 = QPPatchesTest.testCA(model2,1,5) ;
%             ps2 = QPPatchesTest.testPatches(cA2,1,{{5}}) ;
%             dA2 = QPPatchesTest.testDA(cA2,ps2,0) ;
%             compareQP2vsQP3(testCase,dA2)
%         end
        
%         function QP3vsQP2SWIP1x1Heterogeneous(testCase)
%             model2 = QPPatchesTest.testModel(2,[3 3]) ;
%             cA2 = QPPatchesTest.testCA(model2,2,5) ;
%             ps2 = QPPatchesTest.testPatches(cA2,1,{{5}}) ;
%             dA2 = QPPatchesTest.testDA(cA2,ps2,0) ;
%             compareQP2vsQP3(testCase,dA2)
%         end
        
%         function QP3vsQP2SWIP1x2Heterogeneous(testCase)
%             model2 = QPPatchesTest.testModel(2,[4 3]) ;
%             cA2 = QPPatchesTest.testCA(model2,2,[6;7]) ;
%             ps2 = QPPatchesTest.testPatches(cA2,1,{{[6;7]}}) ;
%             dA2 = QPPatchesTest.testDA(cA2,ps2,0) ;
%             compareQP2vsQP3(testCase,dA2)
%         end
        
%         function QP3vsQP2SWIP1x2HeterogeneousSrc(testCase)
%             model2 = QPPatchesTest.testModel(2,[4 3]) ;
%             cA2 = QPPatchesTest.testCA(model2,2,[6;7]) ;
%             ps2 = QPPatchesTest.testPatches(cA2,1,{{[6;7]}}) ;
%             src = QPPatchesTest.testSource(model2) ;
%             dA2 = QPPatchesTest.testDA(cA2,ps2,src) ;
%             compareQP2vsQP3(testCase,dA2)
%         end
        
        function QP3vsQP2ND1x1Homogeneous(testCase)
            model2 = QPPatchesTest.testModel(2,[3 3]) ;
            cA2 = QPPatchesTest.testCA(model2,1,5) ;
            ps2 = QPPatchesTest.testPatches(cA2,2,{{5}}) ;
            dA2 = QPPatchesTest.testDA(cA2,ps2,0) ;
            compareQP2vsQP3(testCase,dA2)
        end
        
        function QP3vsQP2ND1x1Heterogeneous(testCase)
            model2 = QPPatchesTest.testModel(2,[3 3]) ;
            cA2 = QPPatchesTest.testCA(model2,2,5) ;
            ps2 = QPPatchesTest.testPatches(cA2,2,{{5}}) ;
            dA2 = QPPatchesTest.testDA(cA2,ps2,0) ;
            compareQP2vsQP3(testCase,dA2)
        end
        
        function QP3vsQP2ND1x2Heterogeneous(testCase)
            model2 = QPPatchesTest.testModel(2,[4 3]) ;
            cA2 = QPPatchesTest.testCA(model2,2,[6;7]) ;
            ps2 = QPPatchesTest.testPatches(cA2,2,{{[6;7]}}) ;
            dA2 = QPPatchesTest.testDA(cA2,ps2,0) ;
            compareQP2vsQP3(testCase,dA2)
        end
        
        function QP3vsQP2ND1x2HeterogeneousSrc(testCase)
            model2 = QPPatchesTest.testModel(2,[4 3]) ;
            cA2 = QPPatchesTest.testCA(model2,2,[6;7]) ;
            ps2 = QPPatchesTest.testPatches(cA2,2,{{[6;7]}}) ;
            src = QPPatchesTest.testSource(model2) ;
            dA2 = QPPatchesTest.testDA(cA2,ps2,src) ;
            compareQP2vsQP3(testCase,dA2)
        end
        
        function QP3vsQP2Lagrange1x1Homogeneous(testCase)
            model2 = QPPatchesTest.testModel(2,[3 3]) ;
            cA2 = QPPatchesTest.testCA(model2,1,5) ;
            ps2 = QPPatchesTest.testPatches(cA2,3,{{5}}) ;
            dA2 = QPPatchesTest.testDA(cA2,ps2,0) ;
            compareQP2vsQP3(testCase,dA2)
        end
        
        function QP3vsQP2Lagrange1x1Heterogeneous(testCase)
            model2 = QPPatchesTest.testModel(2,[3 3]) ;
            cA2 = QPPatchesTest.testCA(model2,2,5) ;
            ps2 = QPPatchesTest.testPatches(cA2,3,{{5}}) ;
            dA2 = QPPatchesTest.testDA(cA2,ps2,0) ;
            compareQP2vsQP3(testCase,dA2)
        end
        
        function QP3vsQP2Lagrange1x2Heterogeneous(testCase)
            model2 = QPPatchesTest.testModel(2,[4 3]) ;
            cA2 = QPPatchesTest.testCA(model2,2,[6;7]) ;
            ps2 = QPPatchesTest.testPatches(cA2,3,{{[6;7]}}) ;
            dA2 = QPPatchesTest.testDA(cA2,ps2,0) ;
            compareQP2vsQP3(testCase,dA2)
        end
        
        function QP3vsQP2Lagrange1x2HeterogeneousSrc(testCase)
            model2 = QPPatchesTest.testModel(2,[4 3]) ;
            cA2 = QPPatchesTest.testCA(model2,2,[6;7]) ;
            ps2 = QPPatchesTest.testPatches(cA2,3,{{[6;7]}}) ;
            src = QPPatchesTest.testSource(model2) ;
            dA2 = QPPatchesTest.testDA(cA2,ps2,src) ;
            compareQP2vsQP3(testCase,dA2)
        end
        
%         function QP3vsQP2SWIP1x2ND2x1HeterogeneousSrc(testCase)
%             model2 = QPPatchesTest.testModel(2,[6 4]) ;
%             cA2 = QPPatchesTest.testCA(model2,2,[8;9;11;17]) ;
%             ps2 = QPPatchesTest.testPatches(cA2,[1 2],{{[8;9]};{[11;17]}}) ;
%             src = QPPatchesTest.testSource(model2) ;
%             dA2 = QPPatchesTest.testDA(cA2,ps2,src) ;
%             compareQP2vsQP3(testCase,dA2)
%         end
    end
    
    %% Static
    
    methods (Static)
        
        function model = testModel(order,cellNum)
            model = QPModel('order',order,'cellNum',cellNum,'verbose',...
                false,'cellSize',[1 1],'elementSize',[1 1]/10,'tolSVD',1e-6) ;
        end
        
        function cA = testCA(model,type,cells)
            switch type
                case 1 % Homogeneous
                    cA = QPConductivityAssembler.homogeneous(1.5,model) ;
                case 2 % Discs
                    pattern = struct('name',{'uniform','disc'},'value',...
                        {1 9},'size',{[],.2},'center',{[],[.5 .5]},...
                        'offset',{[],[]}) ;
                    cellNb = getCellNb(model) ;
                    dist = {setdiff(1:cellNb,cells)' cells} ;
                    cA = QPConductivityAssembler('model',model,...
                        'patterns',pattern,'patternsTable',[1 1 ; 1 0],...
                        'distribution',dist) ;
            end
            cA = assemble(cA) ;
        end
        
        function src = testSource(model)
            dsz = getDomainSize(model) ;
            src = @(x) exp(-10*((x(:,1)-dsz(1)/2).^2+(x(:,2)-dsz(2)/2).^2)) ;
        end
        
        function ps = testPatches(cA,type,cells)
            cells = formatIndex(getOrder(cA),getCellNum(cA),cells) ;
            Kpatch = QPPatchesTest.subConductivity(getConductivity(cA),cells) ;
            ps = QPPatches('conductivityAssembler',cA,'cells',cells,...
                'useCompression',false,'type',type,'solver',1,...
                'conductivity',Kpatch) ;
        end
        
        function dA = testDA(cA,ps,source) 
            % Typical PBC
            model = getModel(cA) ;
            bc = QPBC('model',model,'type',5,'value',[],'penalty',[],...
                'cells',boundaryCells(model),'isSIP',false,...
                'useAverageWeights',false,'useStabilisationWeights',false) ;
            dA = QPDiffusionAssembler('conductivityAssembler',cA,'bc',bc,...
                'source',source,'penalty',[],'patches',ps,'useCompression',...
                false,'constantNullification','full','useAverageWeights',...
                true,'useStabilisationWeights',true);
            dA = assemble(dA) ;
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
        
        function dAT = QPPatchesTest(tol,errFun)
            if nargin < 2
                errFun = @(a,b) max(-1,norm(full(b-a))/norm(full(a))) ;
                if nargin < 1
                    tol = 1e-10 ;
                end
            end
            dAT.testTol = tol ;
            dAT.testCrit = errFun ;
        end
        
        function [] = compareQP2vsQP3(testCase,dA2)
            ps2 = getPatches(dA2) ;
            model2 = getModel(dA2) ;
            ps2 = assemble(ps2,getMicroOperators(dA2),getPenalty(dA2)) ;
            lhs2 = untensorize(model2,getLHSOperator(ps2),0) ;
            rhs2 = untensorize(model2,getRHSOperator(ps2),0) ;
            model3 = setOrder(model2,3) ;
            cA3 = updateModel(getConductivityAssembler(ps2),model3) ;
            pCells = formatIndex(3,getCellNum(model3),getCells(ps2)) ;
            Kpatch3 = QPPatchesTest.subConductivity(getConductivity(cA3),pCells) ;
            ps3 = setPatchConductivity(ps2,Kpatch3) ;
            ps3 = updateConductivityAssembler(ps3,cA3) ;
            dA3 = setPatches(dA2,ps3) ;
            dA3 = updateModel(dA3,model3) ;
            ps3 = getPatches(dA3) ;
            lhs3 = untensorize(model3,getLHSOperator(ps3),0) ;
            rhs3 = untensorize(model3,getRHSOperator(ps3),0) ;
            errorLHS = testCase.testCrit(lhs2,lhs3) ;
            errorRHS1 = testCase.testCrit(rhs2{1},rhs3{1}) ;
            errorRHS2 = testCase.testCrit(rhs2{2},rhs3{2}) ;
            errorQP = max([errorLHS errorRHS1 errorRHS2]) ;
            testCase.verifyLessThan(errorQP,testCase.testTol) ;
        end
        
    end
end