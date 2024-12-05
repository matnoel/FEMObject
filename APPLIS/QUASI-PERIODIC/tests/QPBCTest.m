classdef QPBCTest < matlab.unittest.TestCase
    
    properties
        iter = 10;
        errFun = @(ref,rel)norm(full(rel-ref))/norm(full(ref)) ;
        tol = 1e-10 ;
    end
    
    methods
        
        function testCase = QPBCTest(iter,errFun,errTol)
            if nargin > 0
                testCase.iter = iter ;
                if nargin > 1
                    testCase.errFun = errFun ;
                    if nargin > 2
                        testCase.tol = errTol ;
                    end
                end
            end
        end
        
        function bool = compareBCOp2FE(testCase,pb)
            [lhsQP,rhsQP] = pb.bc.operators(pb.model,pb.K,true,1) ;
            lhsQP = untensorize(pb.model,lhsQP) ;
            rhsQP = untensorize(pb.model,rhsQP) ;
            
            [~,out] = solveFE(pb) ;
            qpCoord = smooth(pb.model,getDomainCoord(pb.model)) ;
            feCoord = getcoord(out.model.node) ;
            lhsFE = relativeSort(out.operators.lhsBCOp{1},feCoord,qpCoord) ;
            rhsFE = relativeSort(out.operators.rhsBCOp{1},feCoord,qpCoord) ;
            
            
            lError = testCase.errFun(lhsFE,lhsQP) ;
            if pb.bc.type==2
                lError = 0 ; % irrelevant, tested with periodic on 1x1
            end
            if norm(rhsFE)<1e-15
                rError = norm(full(rhsQP)) ;
            else
                rError = testCase.errFun(rhsFE,rhsQP) ;
            end
            
            bool = max(lError,rError)<testCase.tol ;
        end
        
    end
    
    methods (Test)
        
        function construction(testCase)
            generator = {@zeros @ones @rand @randn} ;
            noError = true ;
            for i = 1:testCase.iter
                for j = 1:numel(generator)
                    try
                        QPBC.random([],generator{j}) ;
                    catch
                        noError = false ;
                        break
                    end
                end
            end
            testCase.verifyTrue(noError) ;
        end
        
        function dirichletOperators(testCase)
            allPassed = true ;
            for n = 1:testCase.iter
                model = QPModel.createRandom('verbose',false,...
                    'tolSVD',testCase.tol) ;
                bc = QPBC.random(1) ;
                K = model.scalarField() ;
                pb = QPProblem(model,K,0,bc) ;
                allPassed = allPassed & compareBCOp2FE(testCase,pb) ;
                if ~allPassed
                    break
                end
            end
            testCase.verifyTrue(allPassed)
        end
        
        function neumannOperators(testCase)
            allPassed = true ;
            for n = 1:testCase.iter
                model = QPModel.createRandom('verbose',false,...
                    'tolSVD',testCase.tol) ;
                bc = QPBC.random(2) ;
                K = model.scalarField() ;
                pb = QPProblem(model,K,0,bc) ;
                allPassed = allPassed & compareBCOp2FE(testCase,pb) ;
                if ~allPassed
                    break
                end
            end
            testCase.verifyTrue(allPassed)
        end
        
        function robinOperators(testCase)
            allPassed = true ;
            for n = 1:testCase.iter
                model = QPModel.createRandom('verbose',false,...
                    'tolSVD',testCase.tol) ;
                bc = QPBC.random(3) ;
                K = model.scalarField() ;
                pb = QPProblem(model,K,0,bc) ;
                allPassed = allPassed & compareBCOp2FE(testCase,pb) ;
                if ~allPassed
                    break
                end
            end
            testCase.verifyTrue(allPassed)
        end
        
        function periodicOperators(testCase)
            allPassed = true ;
            for n = 1:testCase.iter
                model = QPModel.createRandom('verbose',false,...
                    'tolSVD',testCase.tol,'cellNum',[1 1]) ;
                bc = QPBC.random(4) ;
                K = model.scalarField() ;
                pb = QPProblem(model,K,0,bc) ;
                allPassed = allPassed & compareBCOp2FE(testCase,pb) ;
                if ~allPassed
                    break
                end
            end
            testCase.verifyTrue(allPassed)
        end
    end
    
    methods (Static)
        
    end
    
end