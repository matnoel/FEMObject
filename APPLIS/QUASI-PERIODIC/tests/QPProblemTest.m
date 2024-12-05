classdef QPProblemTest < matlab.unittest.TestCase
    
    properties
        approxTol = 1e-3 ;
    end

    methods
        
        function testCase = QPProblemTest(tol)
            if nargin > 0
                testCase.approxTol = tol ;
            end
        end
        
        function testApproximation(testCase,pb)
            [~,out] = solve(pb) ;
            % / DEBUG
            disp(out.error)
            % DEBUG \
            testCase.verifyLessThan(out.error,testCase.approxTol)
        end
        
        function compareQP2FE(testCase,pb)
           [ref,outFE] = solveFE(pb) ;
           [rel,~,greedy] = solve(pb) ;
           resError = compareResidual2FE(pb.model,greedy.A*rel,...
               outFE.lhs*ref,getcoord(outFE.model.node)) ;
           % / DEBUG
           errFun = @(ref,rel)norm(full(rel-ref))/norm(full(ref)) ;
           fprintf('Residual error %g - Field error %g\n',resError,...
               errFun(ref,untensorize(pb.model,rel)))
           figure ; plot(pb.model,rel) ; colorbar ;
           figure ; plot(pb.model,ref,'coord',getcoord(outFE.model.node)) ; colorbar ;
           % DEBUG \
           testCase.verifyLessThan(resError,testCase.approxTol) ;
        end
        
    end
    
    methods (Test)
        
        function dirichletPb(testCase)
            model = QPModel.createRandom('verbose',true) ;
            K = QPProblemTest.generateConductivity(model) ;
            src = QPProblemTest.generateSource(model) ;
            bc = QPProblemTest.generateBC(1) ;
            pb = QPProblem(model,K,src,bc,'tolerance',testCase.approxTol) ;
            testCase.testApproximation(pb) ;
        end
        
        function neumannPb(testCase)
            model = QPModel.createRandom('verbose',true) ;
            K = QPProblemTest.generateConductivity(model) ;
            src = QPProblemTest.generateSource(model) ;
            bc = QPProblemTest.generateBC(2) ;
            pb = QPProblem(model,K,src,bc,'tolerance',testCase.approxTol) ;
            testCase.testApproximation(pb) ;            
        end
        
        function robinPb(testCase)
            model = QPModel.createRandom('verbose',true) ;
            K = QPProblemTest.generateConductivity(model) ;
            src = QPProblemTest.generateSource(model) ;
            bc = QPProblemTest.generateBC(3) ;
            pb = QPProblem(model,K,src,bc,'tolerance',testCase.approxTol) ;
            testCase.testApproximation(pb) ;            
        end
        
        function periodicPb(testCase)
            model = QPModel.createRandom('verbose',true,'order',2) ;
            K = QPProblemTest.generateConductivity(model) ;
            src = QPProblemTest.generateSource(model) ;
            bc = QPProblemTest.generateBC(4) ;
            pb = QPProblem(model,K,src,bc,'tolerance',testCase.approxTol) ;
            testCase.testApproximation(pb) ;            
        end
    end
    
    methods (Static)
        
        function K = generateConductivity(model)
            K = model.scalarField(@rand) ;
        end        
        
        function src = generateSource(model)
            src = model.scalarField(@rand) ;
        end
        
        function bc = generateBC(type)
            bc = QPBC.random(type) ;
        end
        
    end 
end