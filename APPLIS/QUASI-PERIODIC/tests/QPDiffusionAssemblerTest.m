classdef QPDiffusionAssemblerTest < matlab.unittest.TestCase
    
    properties
        maxIter
        testTol
        testCrit
    end    
    
    %%
    methods (Test)
        
        function construction(testCase)
            noError = true ;
            for iter = 1:testCase.maxIter
                try
                    QPDiffusionAssembler.createRandom('verbose',false) ;
                catch
                    noError = false ;
                    break
                end
            end
            testCase.verifyTrue(noError) ;
        end
        
        function QP2vsFEMDirichletHomogeneous(testCase)
            model = QPDiffusionAssemblerTest.testModel(2,[1 1]) ;
            cA = QPDiffusionAssemblerTest.testCA(model,1,1.5) ;
            bc = QPDiffusionAssemblerTest.testBC(model,1,2,1) ;
            dA = QPDiffusionAssemblerTest.testDA(cA,bc,0) ;
            compareQP2vsFEM(testCase,dA) ;
        end
        
        function QP2vsFEMDirichletHomogeneousSource(testCase)
            model = QPDiffusionAssemblerTest.testModel(2,[1 1]) ;
            cA = QPDiffusionAssemblerTest.testCA(model,1,1.5) ;
            bc = QPDiffusionAssemblerTest.testBC(model,1,2,1) ;
            src = QPDiffusionAssemblerTest.testSource(model) ;
            dA = QPDiffusionAssemblerTest.testDA(cA,bc,src) ;
            compareQP2vsFEM(testCase,dA) ;
        end
        
        function QP2vsFEMNeumannHomogeneous(testCase)
            model = QPDiffusionAssemblerTest.testModel(2,[1 1]) ;
            cA = QPDiffusionAssemblerTest.testCA(model,1,1.5) ;
            bc = QPDiffusionAssemblerTest.testBC(model,2,2,[]) ;
            dA = QPDiffusionAssemblerTest.testDA(cA,bc,0) ;
            compareQP2vsFEM(testCase,dA) ;
        end
        
        function QP2vsFEMRobinHomogeneous(testCase)
            model = QPDiffusionAssemblerTest.testModel(2,[1 1]) ;
            cA = QPDiffusionAssemblerTest.testCA(model,1,1.5) ;
            bc = QPDiffusionAssemblerTest.testBC(model,3,2,1) ;
            dA = QPDiffusionAssemblerTest.testDA(cA,bc,0) ;
            compareQP2vsFEM(testCase,dA) ;
        end
        
        function QP2vsFEMCauchyHomogeneous(testCase)
            model = QPDiffusionAssemblerTest.testModel(2,[1 1]) ;
            cA = QPDiffusionAssemblerTest.testCA(model,1,1.5) ;
            bc = QPDiffusionAssemblerTest.testBC(model,4,{2 3},1) ;
            dA = QPDiffusionAssemblerTest.testDA(cA,bc,0) ;
            compareQP2vsFEM(testCase,dA) ;
        end
        
%         function QP2vsFEMPeriodicHomogeneous(testCase)
%             model = QPDiffusionAssemblerTest.testModel(2,[1 1]) ;
%             cA = QPDiffusionAssemblerTest.testCA(model,1,1.5) ;
%             bc = QPDiffusionAssemblerTest.testBC(model,5,[],[]) ;
%             dA = QPDiffusionAssemblerTest.testDA(cA,bc,0) ;
%             compareQP2vsFEM(testCase,dA) ;
%         end
        
        function QP3vsQP2PeriodicHomogeneous(testCase)
            model = QPDiffusionAssemblerTest.testModel(2,[3 2]) ;
            cA = QPDiffusionAssemblerTest.testCA(model,1,1.5) ;
            bc = QPDiffusionAssemblerTest.testBC(model,5,[],[]) ;
            dA = QPDiffusionAssemblerTest.testDA(cA,bc,0) ;
            compareQP2vsQP3(testCase,dA) ;
        end
        
        function QP3vsQP2PeriodicHeterogenousNoWeights(testCase)
            model = QPDiffusionAssemblerTest.testModel(2,[3 2]) ;
            cA = QPDiffusionAssemblerTest.testCA(model,2,{1 9}) ;
            bc = QPDiffusionAssemblerTest.testBC(model,5,[],[]) ;
            dA = QPDiffusionAssemblerTest.testDA(cA,bc,0) ;
            dA = setUseStabilisationWeights(dA,false) ;
            dA = setUseAverageWeights(dA,false) ;
            compareQP2vsQP3(testCase,dA) ;
        end
        
        function QP3vsQP2PeriodicHeterogenousNoAverageWeights(testCase)
            model = QPDiffusionAssemblerTest.testModel(2,[3 2]) ;
            cA = QPDiffusionAssemblerTest.testCA(model,2,{1 9}) ;
            bc = QPDiffusionAssemblerTest.testBC(model,5,[],[]) ;
            dA = QPDiffusionAssemblerTest.testDA(cA,bc,0) ;
            dA = setUseStabilisationWeights(dA,true) ;
            dA = setUseAverageWeights(dA,false) ;
            compareQP2vsQP3(testCase,dA) ;
        end
        
        function QP3vsQP2PeriodicHeterogenous(testCase)
            model = QPDiffusionAssemblerTest.testModel(2,[3 2]) ;
            cA = QPDiffusionAssemblerTest.testCA(model,2,{1 9}) ;
            bc = QPDiffusionAssemblerTest.testBC(model,5,[],[]) ;
            dA = QPDiffusionAssemblerTest.testDA(cA,bc,0) ;
            dA = setUseStabilisationWeights(dA,true) ;
            dA = setUseAverageWeights(dA,true) ;
            compareQP2vsQP3(testCase,dA) ;
        end
                
        function QP3vsQP2DirichletHomogeneous(testCase)
            model = QPDiffusionAssemblerTest.testModel(2,[3 2]) ;
            cA = QPDiffusionAssemblerTest.testCA(model,1,1.5) ;
            bc = QPDiffusionAssemblerTest.testBC(model,1,2,1) ;
            dA = QPDiffusionAssemblerTest.testDA(cA,bc,0) ;
            compareQP2vsQP3(testCase,dA) ;
        end
                
        function QP3vsQP2DirichletHomogeneousSource(testCase)
            model = QPDiffusionAssemblerTest.testModel(2,[3 2]) ;
            cA = QPDiffusionAssemblerTest.testCA(model,1,1.5) ;
            bc = QPDiffusionAssemblerTest.testBC(model,1,2,1) ;
            src = QPDiffusionAssemblerTest.testSource(model) ;
            dA = QPDiffusionAssemblerTest.testDA(cA,bc,src) ;
            compareQP2vsQP3(testCase,dA) ;
        end
                
        function QP3vsQP2NeumannHomogeneous(testCase)
            model = QPDiffusionAssemblerTest.testModel(2,[3 2]) ;
            cA = QPDiffusionAssemblerTest.testCA(model,1,1.5) ;
            bc = QPDiffusionAssemblerTest.testBC(model,2,2,[]) ;
            dA = QPDiffusionAssemblerTest.testDA(cA,bc,0) ;
            compareQP2vsQP3(testCase,dA) ;
        end
        
        function QP3vsQP2RobinHomogeneous(testCase)
            model = QPDiffusionAssemblerTest.testModel(2,[3 2]) ;
            cA = QPDiffusionAssemblerTest.testCA(model,1,1.5) ;
            bc = QPDiffusionAssemblerTest.testBC(model,3,2,1) ;
            dA = QPDiffusionAssemblerTest.testDA(cA,bc,0) ;
            compareQP2vsQP3(testCase,dA) ;
        end
        
        function QP3vsQP2CauchyHomogeneous(testCase)
            model = QPDiffusionAssemblerTest.testModel(2,[3 2]) ;
            cA = QPDiffusionAssemblerTest.testCA(model,1,1.5) ;
            bc = QPDiffusionAssemblerTest.testBC(model,4,{2 3},1) ;
            dA = QPDiffusionAssemblerTest.testDA(cA,bc,0) ;
            compareQP2vsQP3(testCase,dA) ;
        end     
         
    end
    
    %% Static
    
    methods (Static)
        
        function model = testModel(order,cellNum)
            model = QPModel('order',order,'cellNum',cellNum,'verbose',...
                false,'cellSize',[1 1],'elementSize',[1 1]/10,'tolSVD',1e-6) ;
        end
        
        function cA = testCA(model,type,value)
            switch type
                case 1 % Homogeneous
                    cA = QPConductivityAssembler.homogeneous(value,model) ;
                case 2 % Discs
                    pattern = struct('name',{'uniform','disc'},'value',...
                        value,'size',{[],.2},'center',{[],[.5 .5]},...
                        'offset',{[],[]}) ;
                    cellNb = getCellNb(model) ;
                    discCell = round(cellNb/2) ;
                    dist = {setdiff(1:cellNb,discCell)' discCell} ;
                    cA = QPConductivityAssembler('model',model,...
                        'patterns',pattern,'patternsTable',[1 1 ; 0 1],...
                        'distribution',dist) ;
            end
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
        
        function dA = testDA(cA,bc,source)
            dA = QPDiffusionAssembler('conductivityAssembler',cA,'bc',bc,...
                'source',source,'penalty',[],'patches',{},'useCompression',...
                false,'constantNullification','full','useAverageWeights',...
                true,'useStabilisationWeights',true);
        end
        
    end
    
    %%
    
    methods
        
        function dAT = QPDiffusionAssemblerTest(tol,errFun,maxIter)
            if nargin < 3
                maxIter = 10 ;
                if nargin < 2
                    errFun = @(a,b) max(-1,norm(full(b-a))/norm(full(a))) ;
                    if nargin < 1
                        tol = 1e-10 ;
                    end
                end
            end
            dAT.testTol = tol ;
            dAT.testCrit = errFun ;
            dAT.maxIter = maxIter ;
        end
        
        function [] = compareQP2vsFEM(testCase,dA)
            dA = assemble(dA) ;
            model = getModel(dA) ;
            qpLHS = untensorize(model,getLHSOperator(dA),0) ;
            qpRHS = untensorize(model,getRHSOperator(dA),0) ;
            pb = QPDiffusionProblem('operatorAssembler',dA) ;
            [~,outFE] = solveFEM(pb) ;
            femCoord = getcoord(getnode(outFE.model)) ;
            qpCoord = smooth(model,getDomainCoord(model)) ;
            femLHS = relativeSort(outFE.lhs,femCoord,qpCoord) ;
            femRHS = relativeSort(outFE.rhs,femCoord,qpCoord) ;
            errorLHS = testCase.testCrit(femLHS,qpLHS) ;
            errorRHS = testCase.testCrit(femRHS,qpRHS) ;
            errorQP2 = max(errorRHS,errorLHS) ;
            testCase.verifyLessThan(errorQP2,testCase.testTol) ;
        end
        
        function [] = compareQP2vsQP3(testCase,dA2)
            dA2 = assemble(dA2) ;
            model2 = getModel(dA2) ;
            lhs2 = untensorize(model2,getLHSOperator(dA2),0) ;
            rhs2 = untensorize(model2,getRHSOperator(dA2),0) ;
            model3 = setOrder(model2,3) ;
            dA3 = updateModel(dA2,model3) ;
            lhs3 = untensorize(model3,getLHSOperator(dA3),0) ;
            rhs3 = untensorize(model3,getRHSOperator(dA3),0) ;
            errorLHS = testCase.testCrit(lhs2,lhs3) ;
            errorRHS = testCase.testCrit(rhs2,rhs3) ;
            errorQP3 = max(errorLHS,errorRHS) ;
            testCase.verifyLessThan(errorQP3,testCase.testTol) ;
        end
        
    end
end