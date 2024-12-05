classdef RackwitzFiesslerSolver
    %
    
    %
    %   A MATLAB class for representing RackwitzFiesslerSolver
    %
    
    properties( SetAccess = public, GetAccess = public )
        maxIterations
        tolerance
        errorCriterion
        display
        displayIterations
    end
    
    methods( Access = public )
        
        % Constructor
        function s = RackwitzFiesslerSolver(varargin)
            % function s = RackwitzFiesslerSolver(varargin)
            % Rackwitz-Fiessler Algorithm
            % s.maxIterations: maximum number of iterations, 100 by default
            % s.tolerance: prescribed tolerance, eps by default
            % s.display: display (true or false), true by default
            % s.displayIterations: display error and stagnation indicators at each step (true or false), false by default
            % s.errorCriterion : error criterion ('designpoint' or 'reliabilityindex'), 'reliabilityindex' by default
            
            expectedErrorCriteria = {'designPoint','reliabilityIndex'};
            
            p = ImprovedInputParser;
            addParameter(p,'maxIterations',100,@isscalar);
            addParameter(p,'tolerance',eps,@isscalar);
            addParameter(p,'errorCriterion','reliabilityIndex',...
                @(x) any(validatestring(x,expectedErrorCriteria)));
            addParameter(p,'display',true,@islogical);
            addParameter(p,'displayIterations',false,@islogical);
            parse(p,varargin{:});
            s = passMatchedArgsToProperties(p,s);
            
        end
        
        function [P,beta,Pf,output] = solve(s,H,gradH,P0,varargin)
            % function [P,beta,Pf,output] = solve(s,H,gradH,P0,varargin)
            % Rackwitz-Fiessler iterative algorithm
            % s: RackwitzFiesslerSolver
            % H: limit state (performance) function
            % gradH: gradient of limit state (performance) function
            % P0: initial design point
            % P: design point
            % beta: reliability index
            % Pf: failure probability
            % output.iter : iteration number
            % output.time : CPU time
            % output.error : error on reliability index beta (if s.errorCriterion = 'reliabilityIndex')
            %                      or design point P (if s.errorCriterion = 'designPoint')
            
            t = tic;
            
            if s.display
                fprintf('\n ---------------------------------------------------------------');
                fprintf('\n ------------ Rackwitz-Fiessler iterative algorithm ------------')
                fprintf('\n ---------------------------------------------------------------\n');
            end
            
            % Initialization
            P = P0;
            beta = norm(P0);
            Pf = normcdf(-beta);
            time = zeros(1,s.maxIterations);
            err = zeros(1,s.maxIterations);
            if s.displayIterations
                fprintf('\nInitialization : P = (%.3f,%.3f)',P(1),P(2));
                fprintf('\n                 beta = %.3f',beta);
                fprintf('\n                 Pf = %.3e\n',Pf);
            end
            
            % Iteration step - Main loop
            for iter=1:s.maxIterations
                
                t_iter = tic;
                
                if nargin(gradH)==2
                    gradHP = gradH(P(1),P(2));
                else
                    gradHP = gradH(P);
                end
                if nargin(H)==2
                    HP = H(P(1),P(2));
                else
                    HP = H(P);
                end

                % Unit normal alpha
                alpha = gradHP./norm(gradHP);
                
                % Reliability index beta
                beta_old = beta;
                % beta = (HP-gradHP'*P)./norm(gradHP);
                beta = HP./norm(gradHP)-alpha'*P;
                
                % Design point P
                P_old = P;
                P = -beta.*alpha;
                
                % Probability failure
                Pf = normcdf(-beta);
                
                switch lower(s.errorCriterion)
                    case 'designpoint'
                        err(iter) = norm(P-P_old);
                    case 'reliabilityindex'
                        err(iter) = abs(beta-beta_old);
                end
                
                time(iter) = toc(t_iter);
                
                % Display
                if s.displayIterations
                    fprintf('\nIteration #%2.d : P = (%.3f,%.3f)',iter,P(1),P(2));
                    fprintf('\n                alpha = (%.3f,%.3f)',alpha(1),alpha(2));
                    fprintf('\n                beta = %.3f',beta);
                    fprintf('\n                Pf = %.3e',Pf);
                    fprintf('\n                error = %.3e',err(iter));
                    fprintf('\n                elapsed time = %f s\n',time(iter));
                end
                % Check for convergence
                if err(iter)<=s.tolerance
                    break
                end
            end
            
            if s.display
                if err(iter)<=s.tolerance
                    fprintf('\nAlgorithm converged at iteration #%d with error = %.3e',iter,err(iter));
                else
                    fprintf('\nAlgorithm stopped at iteration #%d with error = %.3e',iter,err(iter));
                end
                fprintf('\n')
                fprintf('\nElapsed time = %f s\n',toc(t));
            end
            
            % Save outputs
            if nargout>3
                output.iter = iter;
                output.time = time;
                output.error = err;
            end
            
        end

        
    end
    
end
