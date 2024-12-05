classdef IterativeSolver
    %
    
    %
    %   A MATLAB class for representing a multiscale iterative solver
    %
    
    properties
        maxIterations
        tolerance
        relaxation
        updateRelaxationParameter
        errorCriterion
        referenceSolution
        timeSolver
        timeOrder
        display
        displayIterations
    end
    
    methods
        
        % Constructor
        function I = IterativeSolver(varargin)
            % Class IterativeSolver
            %
            % function I = IterativeSolver(varargin)
            % I.maxIterations: maximum number of iterations, 100 by default
            % I.tolerance: prescribed tolerance, eps by default
            % I.relaxation: relaxation (double or 'optimal' or 'optimalApproximation' or 'Aitken'), 1 by default
            % I.updateRelaxationParameter: update relaxation parameter at each iteration (true or false), true by default if rho = 'Aitken', false otherwise
            % I.errorCriterion: error criterion ('none' or 'reference' or 'residual'), 'reference' by default
            % For I.errorCriterion = 'reference',
            % I.referenceSolution: reference solution
            % I.timeSolver: time solver for multiscale problem ([], DGTIMESOLVER, EULERSOLVER, NEWMARKSOLVER), [] by default
            % I.timeOrder: time order for multiscale problem, 0 by default
            % I.display: display error and stagnation indicators at final step (true or false), true by default
            % I.displayIterations: display error and stagnation indicators at each step (true or false), false by default
            
            expectedRelaxationParameters = {'optimal','optimalApproximation','Aitken'};
            expectedErrorCriteria = {'none','reference','residual'};
            
            p = ImprovedInputParser;
            addParameter(p,'maxIterations',100,@isscalar);
            addParameter(p,'tolerance',eps,@isscalar);
            addParameter(p,'relaxation',1,...
                @(x) (isscalar(x) && x>0) || any(validatestring(x,expectedRelaxationParameters)));
            addParameter(p,'updateRelaxationParameter',false,@islogical);
            addParameter(p,'errorCriterion','none',...
                @(x) any(validatestring(x,expectedErrorCriteria)));
            addParameter(p,'referenceSolution',@iscell);
            addParameter(p,'timeSolver',[]);
            addParameter(p,'timeOrder',0,@isscalar);
            addParameter(p,'display',true,@islogical);
            addParameter(p,'displayIterations',false,@islogical);
            parse(p,varargin{:});
            I = passMatchedArgsToProperties(p,I);
            
            if isscalar(I.relaxation)
                I.updateRelaxationParameter = false;
            elseif ischar(I.relaxation) && strcmpi(I.relaxation,'aitken')
                I.updateRelaxationParameter = true;
            end
            
        end
        
        function [U,w,lambda,output,vU,vw,vlambda,aU,aw,alambda] = solve(I,glob,patches,interfaces,varargin)
            % function [U,w,lambda,output] = solve(I,glob,patches,interfaces,varargin)
            % function [U,w,lambda,output,vU,vw,vlambda] = solve(I,glob,patches,interfaces,varargin)
            % function [U,w,lambda,output,vU,vw,vlambda,aU,aw,alambda] = solve(I,glob,patches,interfaces,varargin)
            % Solves multiscale deterministic problem based on overlapping
            % domain decomposition using global-local iterative solver
            %
            % Inputs:
            % I: IterativeSolver
            % glob: Global
            % patches: Patches
            % interfaces: Interfaces
            %
            % Outputs:
            % U: m_U-by-1 doubles containing global solution U
            % w: 1-by-n cell of m_w-by-1 doubles containing local solution w
            % lambda: 1-by-n cell of m_l-by-1 doubles containing Lagrange multiplier lambda
            % output.iteration: 1-by-1 double containing iteration number k
            % output.totalTime: 1-by-1 double containing total CPU time
            % output.time: 1-by-k doubles containing CPU time t
            % output.relaxationParameter: 1-by-k doubles containing relaxation parameter rho
            % output.errorGlobalSolutionInit: 1-by-k doubles containing the error indicator of initial global solution U
            % output.errorLocalSolutionInit: 1-by-n cell of 1-by-k doubles containing the error indicator of initial local solution w
            % output.errorLagrangeMultiplierInit: 1-by-n cell of 1-by-k doubles containing the error indicator of initial Lagrange multiplier lambda
            % output.errorGlobalSolution: 1-by-k doubles containing the error indicator of global solution U
            % output.errorLocalSolution: 1-by-n cell of 1-by-k doubles containing the error indicator of local solution w
            % output.errorLagrangeMultiplier: 1-by-n cell of 1-by-k doubles containing the error indicator of Lagrange multiplier lambda
            % output.stagnationGlobalSolution: 1-by-k doubles containing the stagnation indicator of global solution U
            % output.stagnationLocalSolution: 1-by-n cell of 1-by-k doubles containing the stagnation indicator of local solution w
            % output.stagnationLagrangeMultiplier: 1-by-n cell of 1-by-k doubles containing the stagnation indicator of Lagrange multiplier lambda
            % For first- and second-order time-dependent problems,
            % vU: m_U-by-1 doubles containing global velocity vU
            % vw: 1-by-n cell of m_w-by-1 doubles containing local velocity vw
            % vlambda: 1-by-n cell of m_l-by-1 doubles containing Lagrange multiplier velocity vlambda
            % For second-order time-dependent problems,
            % aU: m_U-by-1 doubles containing global acceleration aU
            % aw: 1-by-n cell of m_w-by-1 doubles containing local acceleration aw
            % alambda: 1-by-n cell of m_l-by-1 doubles containing Lagrange multiplier acceleration alambda
            % where
            % n is the number of patches
            % k is the number of iterations
            % m_U is the dimension of the spatial approximation space of global solution U
            % m_w is the dimension of the spatial approximation space of local solution w
            % m_l is the dimension of the spatial approximation space of Lagrange multiplier lambda
            
            if isscalar(I.relaxation)
                I.updateRelaxationParameter = false;
            elseif ischar(I.relaxation) && strcmpi(I.relaxation,'aitken')
                I.updateRelaxationParameter = true;
            end
            
            if strcmpi(I.errorCriterion,'reference')
                U_ref = I.referenceSolution{1};
                w_ref = I.referenceSolution{2};
                lambda_ref = I.referenceSolution{3};
            else
                [U_ref,w_ref,lambda_ref] = initializeVariables(glob,patches,interfaces);
            end
            
            tTotal = tic;
            
            if I.display || I.displayIterations
                fprintf('\n -------------------------------------------------------');
                fprintf('\n ------------ Global-local iterative solver ------------')
                fprintf('\n -------------------------------------------------------\n');
            end
            
            n = numel(patches);
            
            % Initialization
            [U,w,lambda] = initializeVariables(glob,patches,interfaces);
            if ~isempty(glob.timeSolver)
                if glob.timeOrder>=1
                    [vU,vw,vlambda] = initializeVariables(glob,patches,interfaces);
                end
                if glob.timeOrder>=2
                    [aU,aw,alambda] = initializeVariables(glob,patches,interfaces);
                end
            end
            
            switch lower(I.errorCriterion)
                case 'none'
                    error_U_init = 1;
                    error_w_init = repmat({1},[1,n]);
                    error_lambda_init = repmat({1},[1,n]);
                case 'reference'
                    U_out = glob.P_out*U;
                    error_U_init = norm(U_out - U_ref)/norm(U_ref);
                    error_w_init = cellfun(@(x,xref) norm(x-xref)/norm(xref),w,w_ref,'UniformOutput',false);
                    error_lambda_init = cellfun(@(x,xref) norm(x-xref)/norm(xref),lambda,lambda_ref,'UniformOutput',false);
                case 'residual'
                    error_U_init = 1;
                    error_w_init = repmat({1},[1,n]);
                    error_lambda_init = repmat({1},[1,n]);
            end
            
            % Check for convergence at initialization
            if error_U_init<=I.tolerance
                totalTime = toc(tTotal);
                if I.display || I.displayIterations
                    fprintf('\nAlgorithm converged at initialization with error = %.3e w.r.t. U',error_U_init);
                    for k=1:n
                        fprintf('\n                                          error = %.3e w.r.t. w{%u}',error_w_init{k},k);
                        fprintf('\n                                          error = %.3e w.r.t. lambda{%u}',error_lambda_init{k},k);
                    end
                    fprintf('\n')
                    fprintf('\nElapsed time = %f s\n',totalTime);
                end
                output.iteration = 0;
                output.totalTime = totalTime;
                output.errorGlobalSolutionInit = error_U_init;
                output.errorLocalSolutionInit = error_w_init;
                output.errorLagrangeMultiplierInit = error_lambda_init;
                return
            elseif I.displayIterations
                fprintf('\nInitialization : error = %.3e w.r.t. U',error_U_init);
                for k=1:n
                    fprintf('\n                 error = %.3e w.r.t. w{%u}',error_w_init{k},k);
                    fprintf('\n                 error = %.3e w.r.t. lambda{%u}',error_lambda_init{k},k);
                end
                fprintf('\n');
            end
            
            [error_U,error_w,error_lambda]...
                = initializeIndicators(I.maxIterations,patches);
            [stagnation_U,stagnation_w,stagnation_lambda]...
                = initializeIndicators(I.maxIterations,patches);
            time = initializeIndicators(I.maxIterations,patches);
            relaxationParameter = initializeIndicators(I.maxIterations,patches);
            
            V = U;
            U_old = U;
            if ~isempty(glob.timeSolver)
                if glob.timeOrder>=1
                    vV = vU;
                end
                if glob.timeOrder>=2
                    aV = aU;
                end
            end
            
            % Relaxation parameter rho
            if ~I.updateRelaxationParameter
                rho = I.calcRelaxationParameter(glob,patches,interfaces,V,w);
            else
                rho = 1;
            end
            
            % Iteration step - Main loop
            for iter=1:I.maxIterations
                
                tIter = tic;
                
                % Global problem
                % Global solution V
                V_old = V;
                if isempty(glob.timeSolver)
                    % Time-independent problem
                    V = glob.solve(interfaces,lambda,U,V);
                else
                    % Time-dependent problem
                    if glob.timeOrder==1
                        [V,~,vV] = glob.solve(interfaces,lambda,U,V,vU,vV);
                    elseif glob.timeOrder==2
                        [V,~,vV,aV] = glob.solve(interfaces,lambda,U,V,vU,vV,aU,aV);
                    end
                end
                
                % Relaxation step
                % Relaxation parameter rho
                if I.updateRelaxationParameter
                    rho = I.calcRelaxationParameter(glob,patches,interfaces,V,w,V_old,U,U_old,rho,iter);
                end
                relaxationParameter(iter) = rho;
                
                % Global iterate U
                U_old = U;
                U = rho*V + (1-rho)*U;
                if ~isempty(glob.timeSolver)
                    if glob.timeOrder>=1
                        vU = rho*vV + (1-rho)*vU;
                    end
                    if glob.timeOrder>=2
                        aU = rho*aV + (1-rho)*aU;
                    end
                end
                
                stagnation_U(iter) = norm(U-U_old)/norm(U);
                
                switch lower(I.errorCriterion)
                    case 'none'
                        error_U(iter) = 1;
                    case 'reference'
                        U_out = glob.P_out*U;
                        error_U(iter) = norm(U_out-U_ref)/norm(U_ref);
                    case 'residual'
                        error_U(iter) = 1;
                end
                
                % Local problems
                % Local solutions (w{k},lambda{k}) without change of variable
                %                             z{k} with change of variable
                w_old = w;
                lambda_old = lambda;
                patch = patches.patches;
                interface = interfaces.interfaces;
                for k=1:n
                    if isempty(patch{k}.timeSolver)
                        [w{k},lambda{k}] = patch{k}.solve(interface{k},U,w{k},lambda{k});
                    else
                        if patch{k}.timeOrder==1
                            [w{k},lambda{k},~,vw{k},vlambda{k}] = patch{k}.solve(interface{k},U,w{k},lambda{k},vU,vw{k},vlambda{k});
                        elseif patch{k}.timeOrder==2
                            [w{k},lambda{k},~,vw{k},vlambda{k},aw{k},alambda{k}] = patch{k}.solve(interface{k},U,w{k},lambda{k},vU,vw{k},vlambda{k},aU,aw{k},alambda{k});
                        end
                    end
                    
                    stagnation_w{k}(iter) = norm(w{k}-w_old{k})/norm(w{k});
                    stagnation_lambda{k}(iter) = norm(lambda{k}-lambda_old{k})/norm(lambda{k});
                    
                    switch lower(I.errorCriterion)
                        case 'none'
                            error_w{k}(iter) = 1;
                            error_lambda{k}(iter) = 1;
                        case 'reference'
                            error_w{k}(iter) = norm(w{k}-w_ref{k})/norm(w_ref{k});
                            error_lambda{k}(iter) = norm(lambda{k}-lambda_ref{k})/norm(lambda_ref{k});
                        case 'residual'
                            error_w{k}(iter) = 1;
                            error_w{k}(iter) = 1;
                    end
                end
                
                time(iter) = toc(tIter);
                
                % Check for convergence
                if error_U(iter)<=I.tolerance
                    break
                else
                    if I.displayIterations
                        fprintf('\nIteration #%2.d : stagnation = %.3e, error = %.3e w.r.t. U',iter,stagnation_U(iter),error_U(iter));
                        for k=1:n
                            fprintf('\n                stagnation = %.3e, error = %.3e w.r.t. w{%u}',stagnation_w{k}(iter),error_w{k}(iter),k);
                            fprintf('\n                stagnation = %.3e, error = %.3e w.r.t. lambda{%u}',stagnation_lambda{k}(iter),error_lambda{k}(iter),k);
                        end
                        fprintf('\n                elapsed time = %f s\n',time(iter));
                    end
                end
            end
            
            totalTime = toc(tTotal);
            
            if I.display || I.displayIterations
                if error_U(iter)<=I.tolerance
                    str = sprintf('Algorithm converged at iteration #%d with ',iter);
                    fprintf(['\n' str 'error = %.3e w.r.t. U'],error_U(iter));
                else
                    str = sprintf('Algorithm stopped at iteration #%d with ',iter);
                    fprintf(['\n' str 'error = %.3e w.r.t. U'],error_U(iter));
                end
                for k=1:n
                    fprintf(['\n' repmat(sprintf(' '),1,length(str)) 'error = %.3e w.r.t. w{%u}'],error_w{k}(iter),k);
                    fprintf(['\n' repmat(sprintf(' '),1,length(str)) 'error = %.3e w.r.t. lambda{%u}'],error_lambda{k}(iter),k);
                end
                fprintf('\n')
                fprintf('\nElapsed time = %f s\n',totalTime);
            end
            
            % Save outputs
            output.iteration = iter;
            output.totalTime = totalTime;
            output.time = time;
            output.relaxationParameter = relaxationParameter;
            
            output.errorGlobalSolutionInit = error_U_init;
            output.errorLocalSolutionInit = error_w_init;
            output.errorLagrangeMultiplierInit = error_lambda_init;
            
            output.errorGlobalSolution = error_U;
            output.errorLocalSolution = error_w;
            output.errorLagrangeMultiplier = error_lambda;
            
            output.stagnationGlobalSolution = stagnation_U;
            output.stagnationLocalSolution = stagnation_w;
            output.stagnationLagrangeMultiplier = stagnation_lambda;
            
        end
        
%         function [U,w,lambda,output,vU,vw,vlambda] = dsolve(I,glob,patches,interfaces,varargin)
%             % function [U,w,lambda,output,vU,vw,vlambda] = dsolve(I,glob,patches,interfaces,varargin)
%             % Solves multiscale deterministic problem based on overlapping
%             % domain decomposition using global-local iterative solver
%             %
%             % Inputs:
%             % I: IterativeSolver
%             % glob: Global
%             % patches: Patches
%             % interfaces: Interfaces
%             %
%             % Outputs:
%             % U: m_U-by-1 doubles containing global solution U
%             % w: 1-by-n cell of m_w-by-1 doubles containing local solution w
%             % lambda: 1-by-n cell of m_l-by-1 doubles containing Lagrange multiplier lambda
%             % output.iteration: 1-by-1 double containing iteration number k
%             % output.totalTime: 1-by-1 double containing total CPU time
%             % output.time: 1-by-k doubles containing CPU time t
%             % output.relaxationParameter: 1-by-k doubles containing relaxation parameter rho
%             % output.errorGlobalSolutionInit: 1-by-k doubles containing the error indicator of initial global solution U
%             % output.errorLocalSolutionInit: 1-by-n cell of 1-by-k doubles containing the error indicator of initial local solution w
%             % output.errorLagrangeMultiplierInit: 1-by-n cell of 1-by-k doubles containing the error indicator of initial Lagrange multiplier lambda
%             % output.errorGlobalSolution: 1-by-k doubles containing the error indicator of global solution U
%             % output.errorLocalSolution: 1-by-n cell of 1-by-k doubles containing the error indicator of local solution w
%             % output.errorLagrangeMultiplier: 1-by-n cell of 1-by-k doubles containing the error indicator of Lagrange multiplier lambda
%             % output.stagnationGlobalSolution: 1-by-k doubles containing the stagnation indicator of global solution U
%             % output.stagnationLocalSolution: 1-by-n cell of 1-by-k doubles containing the stagnation indicator of local solution w
%             % output.stagnationLagrangeMultiplier: 1-by-n cell of 1-by-k doubles containing the stagnation indicator of Lagrange multiplier lambda
%             % vU: m_U-by-1 doubles containing global velocity vU
%             % vw: 1-by-n cell of m_w-by-1 doubles containing local velocity vw
%             % vlambda: 1-by-n cell of m_l-by-1 doubles containing Lagrange multiplier velocity vlambda
%             % where
%             % n is the number of patches
%             % k is the number of iterations
%             % m_U is the dimension of the spatial approximation space of global solution U
%             % m_w is the dimension of the spatial approximation space of local solution w
%             % m_l is the dimension of the spatial approximation space of Lagrange multiplier lambda
%             
%             if isscalar(I.relaxation)
%                 I.updateRelaxationParameter = false;
%             elseif ischar(I.relaxation) && strcmpi(I.relaxation,'aitken')
%                 I.updateRelaxationParameter = true;
%             end
%             
%             if strcmpi(I.errorCriterion,'reference')
%                 Ut_ref = I.referenceSolution{1};
%                 wt_ref = I.referenceSolution{2};
%                 lambdat_ref = I.referenceSolution{3};
%             end
%             
%             tTotal = tic;
%             
%             display_ = getparam(I.timeSolver,'display');
%             if isa(I.timeSolver,'EULERTIMESOLVER')
%                 eulertype = getparam(I.timeSolver,'eulertype');
%             end
%             if display_
%                 if isa(I.timeSolver,'EULERTIMESOLVER')
%                     fprintf('\n -------------------------------------------');
%                     fprintf('\n ------------ Euler time solver ------------')
%                     fprintf('\n -------------------------------------------\n');
%                 elseif isa(I.timeSolver,'DGTIMESOLVER')
%                     fprintf('\n ---------------------------------------');
%                     fprintf('\n ------------ DGtime solver ------------')
%                     fprintf('\n ---------------------------------------\n');
%                 end
%             end
%             
%             n = numel(patches);
%             [glob,patches] = initializeRightHandSide(I,glob,patches);
%             
%             patch = patches.patches;
%             interface = interfaces.interfaces;
%             
%             T = gettimemodel(I.timeSolver);
%             % t = gett(T);
%             nt = getnt(T);
%             dt = getdt(T);
%             p  = getp(I.timeSolver);
%             
%             % Initialization
%             Ut = cell(1,length(T));
%             vUt = cell(1,length(T));
%             Ut{1} = initializeSolution(glob);
%             if isa(I.timeSolver,'EULERTIMESOLVER') && strcmpi(eulertype,'implicit')
%                 sz_U = getnbddlfree(glob.S);
%                 vUt{1} = zeros(sz_U,1);
%             end
%             
%             wt = cell(1,n);
%             lambdat = cell(1,n);
%             vwt = cell(1,n);
%             vlambdat = cell(1,n);
%             for k=1:n
%                 wt{k} = cell(1,length(T));
%                 lambdat{k} = cell(1,length(T));
%                 vwt{k} = cell(1,length(T));
%                 vlambdat{k} = cell(1,length(T));
%                 wt{k}{1} = initializeSolution(patch{k});
%                 lambdat{k}{1} = initializeSolution(interface{k});
%                 if isa(I.timeSolver,'EULERTIMESOLVER') && strcmpi(eulertype,'implicit')
%                     sz_w = getnbddl(patch{k}.S);
%                     sz_lambda = getnbddl(interface{k}.S);
%                     vwt{k}{1} = zeros(sz_w,1);
%                     vlambdat{k}{1} = zeros(sz_lambda,1);
%                 end
%             end
%             
%             outputt = cell(1,length(T));
%             
%             if display_
%                 if isa(I.timeSolver,'EULERTIMESOLVER')
%                     fprintf('Euler %s : ',eulertype);
%                 elseif isa(I.timeSolver,'DGTIMESOLVER')
%                     fprintf('DG order %d : ',p);
%                 end
%             end
%             for i=1:nt
%                 if display_
%                     pourcentage(i,nt)
%                 end
%                 
%                 if strcmpi(I.errorCriterion,'reference')
%                     U_ref = Ut_ref{i+1}; % U_ref = getmatrixatstep(Ut_ref,i+1);
%                     w_ref = cellfun(@(patch) wt_ref{patch.number}{i+1},patch,'UniformOutput',false); % w_ref = cellfun(@(patch) getmatrixatstep(wt_ref{patch.number},i+1),patch,'UniformOutput',false);
%                     lambda_ref = cellfun(@(interface) lambdat_ref{interface.number}{i+1},interface,'UniformOutput',false); % lambda_ref = cellfun(@(interface) getmatrixatstep(lambdat_ref{interface.number},i+1),interface,'UniformOutput',false);
%                 else
%                     [U_ref,w_ref,lambda_ref] = initializeVariables(glob,patches,interfaces);
%                 end
%                 
%                 % Initialization
%                 [U,w,lambda] = initializeVariables(glob,patches,interfaces);
%                 [vU,vw,vlambda] = initializeVariables(glob,patches,interfaces);
%                 
%                 switch lower(I.errorCriterion)
%                     case 'none'
%                         error_U_init = 1;
%                         error_w_init = repmat({1},[1,n]);
%                         error_lambda_init = repmat({1},[1,n]);
%                     case 'reference'
%                         U_out = glob.P_out*U;
%                         error_U_init = norm(U_out - U_ref)/norm(U_ref);
%                         error_w_init = cellfun(@(x,xref) norm(x-xref)/norm(xref),w,w_ref,'UniformOutput',false);
%                         error_lambda_init = cellfun(@(x,xref) norm(x-xref)/norm(xref),lambda,lambda_ref,'UniformOutput',false);
%                     case 'residual'
%                         error_U_init = 1;
%                         error_w_init = repmat({1},[1,n]);
%                         error_lambda_init = repmat({1},[1,n]);
%                 end
%                 
%                 % Check for convergence at initialization
%                 if error_U_init<=I.tolerance
%                     if I.display || I.displayIterations
%                         fprintf('\nAlgorithm converged at initialization with error = %.3e w.r.t. U',error_U_init);
%                         for k=1:n
%                             fprintf('\n                                          error = %.3e w.r.t. w{%u}',error_w_init{k},k);
%                             fprintf('\n                                          error = %.3e w.r.t. lambda{%u}',error_lambda_init{k},k);
%                         end
%                         fprintf('\n')
%                         fprintf('\nElapsed time = %f s\n',toc(tTotal));
%                     end
%                     outputt{i}.iteration = 0;
%                     outputt{i}.errorGlobalSolutionInit = error_U_init;
%                     outputt{i}.errorLocalSolutionInit = error_w_init;
%                     outputt{i}.errorLagrangeMultiplierInit = error_lambda_init;
%                     return
%                 elseif I.displayIterations
%                     fprintf('\nInitialization : error = %.3e w.r.t. U',error_U_init);
%                     for k=1:n
%                         fprintf('\n                 error = %.3e w.r.t. w{%u}',error_w_init{k},k);
%                         fprintf('\n                 error = %.3e w.r.t. lambda{%u}',error_lambda_init{k},k);
%                     end
%                     fprintf('\n');
%                 end
%                 
%                 [error_U,error_w,error_lambda]...
%                     = initializeIndicators(I.maxIterations,patches);
%                 [stagnation_U,stagnation_w,stagnation_lambda]...
%                     = initializeIndicators(I.maxIterations,patches);
%                 time = initializeIndicators(I.maxIterations,patches);
%                 relaxationParameter = initializeIndicators(I.maxIterations,patches);
%                 
%                 V = U;
%                 U_old = U;
%                 vV = vU;
%                 
%                 % Relaxation parameter rho
%                 if ~I.updateRelaxationParameter
%                     rho = I.calcRelaxationParameter(glob,patches,interfaces,V,w);
%                 else
%                     rho = 1;
%                 end
%                 
%                 % Iteration step - Main loop
%                 for iter=1:I.maxIterations
%                     
%                     tIter = tic;
%                     
%                     % Global problem
%                     % Global solution V
%                     V_old = V;
%                     %[V,~,vV] = glob.dsolve(interfaces,dt(i),Vt{i},lambda,U,V,vU,vV);
%                     
%                     % Global soluton V
%                     if isa(glob.timeSolver,'EULERTIMESOLVER')
%                         switch lower(eulertype)
%                             case 'explicit'
%                                 b = glob.b_out{i};
%                                 for k=1:n
%                                     B_glob = interface{k}.P_glob'*interface{k}.M;
%                                     b = b - B_glob*lambdat{k}{i} + glob.A_in{k}*Ut{i} + glob.M_in{k}*vUt{i};
%                                 end
%                                 if glob.increment && (isempty(glob.solver) || ~isanlsolver(glob.solver))
%                                     b = b - glob.A*V_old - glob.M*vV_old;
%                                 end
%                                 M = glob.M/dt(i);
%                                 bi = b + glob.M*Vt{i}/dt(i) - glob.A*Vt{i};
%                                 [V,~] = solve(M,bi);
%                                 vV = (V - Vt{i})/dt(i);
%                                 % vV = solve(glob.M,b-glob.A*Vt{i});
%                             case 'implicit'
%                                 b = glob.b_out{i+1};
%                                 for k=1:n
%                                     B_glob = interface{k}.P_glob'*interface{k}.M;
%                                     b = b - B_glob*lambda_old{k} + glob.A_in{k}*U_old + glob.M_in{k}*vU_old;
%                                 end
%                                 if glob.increment && (isempty(glob.solver) || ~isanlsolver(glob.solver))
%                                     b = b - glob.A*V_old - glob.M*vV_old;
%                                 end
%                                 M = glob.M/dt(i)+glob.A;
%                                 bi = b + glob.M*Vt{i}/dt(i);
%                                 [V,~] = solve(M,bi);
%                                 vV = (V - Vt{i})/dt(i);
%                                 % vV = solve(glob.M,b-glob.A*V);
%                         end
%                     elseif isa(glob.timeSolver,'DGTIMESOLVER')
%                         
%                     end
%                     
%                     if glob.increment && (isempty(glob.solver) || ~isanlsolver(glob.solver))
%                         V = V_old + V;
%                         vV = vV_old + vV;
%                     end
%                     
%                     % Relaxation step
%                     % Relaxation parameter rho
%                     if I.updateRelaxationParameter
%                         rho = I.calcRelaxationParameter(glob,patches,interfaces,V,w,V_old,U,U_old,rho,iter);
%                     end
%                     relaxationParameter(iter) = rho;
%                     
%                     % Global iterate U
%                     U_old = U;
%                     U = rho*V + (1-rho)*U;
%                     vU = rho*vV + (1-rho)*vU;
%                     
%                     stagnation_U(iter) = norm(U-U_old)/norm(U);
%                     
%                     switch lower(I.errorCriterion)
%                         case 'none'
%                             error_U(iter) = 1;
%                         case 'reference'
%                             U_out = glob.P_out*U;
%                             error_U(iter) = norm(U_out-U_ref)/norm(U_ref);
%                         case 'residual'
%                             error_U(iter) = 1;
%                     end
%                     
%                     % Local problems
%                     % Local solutions (w{k},lambda{k}) without change of variable
%                     %                             z{k} with change of variable
%                     w_old = w;
%                     lambda_old = lambda;
%                     patch = patches.patches;
%                     interface = interfaces.interfaces;
%                     for k=1:n
%                         % [w{k},lambda{k},~,vw{k},vlambda{k}] = patch{k}.dsolve(interface{k},U,w{k},lambda{k},vU,vw{k},vlambda{k});
%                         
%                         
%                         
%                         stagnation_w{k}(iter) = norm(w{k}-w_old{k})/norm(w{k});
%                         stagnation_lambda{k}(iter) = norm(lambda{k}-lambda_old{k})/norm(lambda{k});
%                         
%                         switch lower(I.errorCriterion)
%                             case 'none'
%                                 error_w{k}(iter) = 1;
%                                 error_lambda{k}(iter) = 1;
%                             case 'reference'
%                                 error_w{k}(iter) = norm(w{k}-w_ref{k})/norm(w_ref{k});
%                                 error_lambda{k}(iter) = norm(lambda{k}-lambda_ref{k})/norm(lambda_ref{k});
%                             case 'residual'
%                                 error_w{k}(iter) = 1;
%                                 error_w{k}(iter) = 1;
%                         end
%                     end
%                     
%                     time(iter) = toc(tIter);
%                     
%                     % Check for convergence
%                     if error_U(iter)<=I.tolerance
%                         break
%                     else
%                         if I.displayIterations
%                             fprintf('\nIteration #%2.d : stagnation = %.3e, error = %.3e w.r.t. U',iter,stagnation_U(iter),error_U(iter));
%                             for k=1:n
%                                 fprintf('\n                stagnation = %.3e, error = %.3e w.r.t. w{%u}',stagnation_w{k}(iter),error_w{k}(iter),k);
%                                 fprintf('\n                stagnation = %.3e, error = %.3e w.r.t. lambda{%u}',stagnation_lambda{k}(iter),error_lambda{k}(iter),k);
%                             end
%                             fprintf('\n                elapsed time = %f s\n',time(iter));
%                         end
%                     end
%                 end
%                 
%                 totalTime = toc(tTotal);
%                 
%                 if I.display || I.displayIterations
%                     if error_U(iter)<=I.tolerance
%                         str = sprintf('Algorithm converged at iteration #%d with ',iter);
%                         fprintf(['\n' str 'error = %.3e w.r.t. U'],error_U(iter));
%                     else
%                         str = sprintf('Algorithm stopped at iteration #%d with ',iter);
%                         fprintf(['\n' str 'error = %.3e w.r.t. U'],error_U(iter));
%                     end
%                     for k=1:n
%                         fprintf(['\n' repmat(sprintf(' '),1,length(str)) 'error = %.3e w.r.t. w{%u}'],error_w{k}(iter),k);
%                         fprintf(['\n' repmat(sprintf(' '),1,length(str)) 'error = %.3e w.r.t. lambda{%u}'],error_lambda{k}(iter),k);
%                     end
%                     fprintf('\n')
%                     fprintf('\nElapsed time = %f s\n',totalTime);
%                 end
%                 
%                 Ut{i+1} = U;
%                 for k=1:n
%                     wt{k}{i+1} = w;
%                     lambdat{k}{i+1} = lambda;
%                 end
%                 switch eulertype
%                     case 'explicit'
%                         vUt{i} = vU;
%                         for k=1:n
%                             vwt{k}{i} = vw;
%                             vlambdat{k}{i} = vlambda;
%                         end
%                     case 'implicit'
%                         vUt{i+1} = vU;
%                         for k=1:n
%                             vwt{k}{i+1} = vw;
%                             vlambdat{k}{i+1} = vlambda;
%                         end
%                 end
%                 
%                 % Save outputs
%                 outputt{i}.iteration = iter;
%                 outputt{i}.totalTime = totalTime;
%                 outputt{i}.time = time;
%                 outputt{i}.relaxationParameter = relaxationParameter;
%                 
%                 outputt{i}.errorGlobalSolutionInit = error_U_init;
%                 outputt{i}.errorLocalSolutionInit = error_w_init;
%                 outputt{i}.errorLagrangeMultiplierInit = error_lambda_init;
%                 
%                 outputt{i}.errorGlobalSolution = error_U;
%                 outputt{i}.errorLocalSolution = error_w;
%                 outputt{i}.errorLagrangeMultiplier = error_lambda;
%                 
%                 outputt{i}.stagnationGlobalSolution = stagnation_U;
%                 outputt{i}.stagnationLocalSolution = stagnation_w;
%                 outputt{i}.stagnationLagrangeMultiplier = stagnation_lambda;
%             end
%             
%             if isa(I.timeSolver,'EULERTIMESOLVER') && strcmpi(eulertype,'implicit')
%                 sz_U = getnbddlfree(glob.S);
%                 vUt{nt+1} = zeros(sz_U,1);
%                 for k=1:n
%                     sz_w = getnbddl(patch{k}.S);
%                     sz_lambda = getnbddl(interface{k}.S);
%                     vwt{k}{nt+1} = zeros(sz_w,1);
%                     vlambdat{k}{nt+1} = zeros(sz_lambda,1);
%                 end
%             end
%             
%             totalTime = toc(tTotal);
%             if display_
%                 fprintf('Elapsed time = %f\n',totalTime)
%             end
%             
%             try
%                 Ut = horzcat(Ut{:});
%                 vUt = horzcat(vUt{:});
%             catch
%                 warning('cell array')
%             end
%             Ut = TIMEMATRIX(Ut,T,[sz_U,1]);
%             vUt = TIMEMATRIX(vUt,T,[sz_U,1]);
%             
%         end
        
        function [U,w,lambda,output] = solveRandom(I,glob,patches,interfaces,s,bases,ls,rv)
            % function [U,w,lambda,output] = solveRandom(I,glob,patches,interfaces,s,bases,ls,rv)
            % Solves multiscale stochastic problem based on overlapping
            % domain decomposition using global-local iterative solver
            %
            % Inputs:
            % I: IterativeSolver
            % glob: Global
            % patches: Patches
            % interfaces: Interfaces
            % s: AdaptiveSparseTensorAlgorithm
            % bases: FunctionalBases
            % ls: LinearModelLearningSquareLoss
            % rv: RandomVector or RandomVariable (optional)
            %
            % Outputs:
            % U: FunctionalBasisArray of size m_U-by-p of global solution U
            % w: 1-by-n cell of FunctionalBasisArray of size m_w-by-p of local solution w
            % lambda: 1-by-n cell of FunctionalBasisArray of size m_l-by-p of Lagrange multiplier lambda
            % output.iteration: 1-by-1 double containing iteration number k
            % output.totalTime: 1-by-1 double containing total CPU time
            % output.time: 1-by-k doubles containing CPU time t
            % output.relaxationParameter: 1-by-k doubles containing relaxation parameter rho
            % output.errorGlobalSolutionInit: 1-by-k doubles containing the error indicator of initial global solution U
            % output.errorLocalSolutionInit: 1-by-n cell of 1-by-k doubles containing the error indicator of initial local solution w
            % output.errorLagrangeMultiplierInit: 1-by-n cell of 1-by-k doubles containing the error indicator of initial Lagrange multiplier lambda
            % output.errorGlobalSolution: 1-by-k doubles containing the error indicator of global solution U
            % output.errorLocalSolution: 1-by-n cell of 1-by-k doubles containing the error indicator of local solution w
            % output.errorLagrangeMultiplier: 1-by-n cell of 1-by-k doubles containing the error indicator of Lagrange multiplier lambda
            % output.stagnationGlobalSolution: 1-by-k doubles containing the stagnation indicator of global solution U
            % output.stagnationLocalSolution: 1-by-n cell of 1-by-k doubles containing the stagnation indicator of local solution w
            % output.stagnationLagrangeMultiplier: 1-by-n cell of 1-by-k doubles containing the stagnation indicator of Lagrange multiplier lambda
            % output.CVErrorLocalSolution: 1-by-n cell of 1-by-k cell of 1-by-m_w doubles containing the cross-validation error of local soluton w
            % output.CVErrorLagrangeMultiplier: 1-by-n cell of 1-by-k cell of 1-by-m_l doubles containing the cross-validation error of Lagrange multiplier lambda
            % output.dimBasisLocalSolution: 1-by-n cell of 1-by-k doubles containing the dimension of the stochastic approximation space of local solution w
            % output.dimBasisLagrangeMultiplier: 1-by-n cell of 1-by-k doubles containing the dimension of the stochastic approximation space of Lagrange multiplier lambda
            % output.nbSamples: 1-by-n cell of 1-by-k doubles containing sample size N
            % output.SamplesLocalSolution: 1-by-n cell of N-by-m_w-by-p doubles containing the evaluations of local solution w
            % output.SamplesLagrangeMultiplier: 1-by-n cell of N-by-m_l-by-p doubles containing the evaluations of Lagrange multiplier lambda
            % where
            % n is the number of patches
            % k is the number of iterations
            % N is the number of samples
            % m_U is the dimension of the spatial approximation space of global solution U
            % m_w is the dimension of the spatial approximation space of local solution w
            % m_l is the dimension of the spatial approximation space of Lagrange multiplier lambda
            % p is the dimension of the time approximation space of multiscale problem
            
            d = length(bases);
            rvb = getRandomVector(bases);
            if nargin<8 || isempty(rv)
                rv = rvb;
            elseif isa(rv,'RandomVariable')
                rv = RandomVector(rv,d);
            end
            
            if isscalar(I.relaxation)
                I.updateRelaxationParameter = false;
            elseif ischar(I.relaxation) && strcmpi(I.relaxation,'aitken')
                I.updateRelaxationParameter = true;
            end
            
            if strcmpi(I.errorCriterion,'reference')
                U_ref = I.referenceSolution{1};
                w_ref = I.referenceSolution{2};
                lambda_ref = I.referenceSolution{3};
            else
                [U_ref,w_ref,lambda_ref] = initializeRandomVariables(glob,patches,interfaces,bases);
            end
            
            tTotal = tic;
            
            if I.display || I.displayIterations
                fprintf('\n -------------------------------------------------------');
                fprintf('\n ------------ Global-local iterative solver ------------')
                fprintf('\n -------------------------------------------------------\n');
            end
            
            n = numel(patches);
            
            % Initialization
            [U,w,lambda] = initializeRandomVariables(glob,patches,interfaces,bases);
            if ~isempty(glob.timeSolver)
                if glob.timeOrder>=1
                    [vU,vw,vlambda] = initializeRandomVariables(glob,patches,interfaces,bases);
                end
                if glob.timeOrder>=2
                    [aU,aw,alambda] = initializeRandomVariables(glob,patches,interfaces,bases);
                end
            end
            
            switch lower(I.errorCriterion)
                case 'none'
                    error_U_init = 1;
                    error_w_init = repmat({1},[1,n]);
                    error_lambda_init = repmat({1},[1,n]);
                case 'reference'
                    indices = addIndices(U.basis.indices,U_ref.basis.indices);
                    basis = SparseTensorProductFunctionalBasis(bases,indices);
                    U_proj = U.projection(basis);
                    U_ref_proj = U_ref.projection(basis);
                    if isempty(glob.timeSolver)
                        % Time-independent problem
                        U_out_proj = U_proj*glob.P_out';
                    else
                        % Time-dependent problem
                        U_out_proj = U_proj;
                        data = permute(U_proj.data,[2 1 3]);
                        sz = size(data);
                        data = glob.P_out*data(:,:);
                        data = reshape(data,[size(data,1),sz(2:3)]);
                        U_out_proj.data = permute(data,[2 1 3]);
                    end
                    error_U_init = norm(U_out_proj - U_ref_proj)/norm(U_ref_proj);
                    
                    indices = cellfun(@(x,y) addIndices(x.basis.indices,y.basis.indices),w,w_ref,'UniformOutput',false);
                    basis = cellfun(@(x) SparseTensorProductFunctionalBasis(bases,x),indices,'UniformOutput',false);
                    w_proj = cellfun(@(x,H) x.projection(H),w,basis,'UniformOutput',false);
                    w_ref_proj = cellfun(@(x,H) x.projection(H),w_ref,basis,'UniformOutput',false);
                    error_w_init = cellfun(@(x,xref) norm(x-xref)/norm(xref),w_proj,w_ref_proj,'UniformOutput',false);
                    
                    indices = cellfun(@(x,y) addIndices(x.basis.indices,y.basis.indices),lambda,lambda_ref,'UniformOutput',false);
                    basis = cellfun(@(x) SparseTensorProductFunctionalBasis(bases,x),indices,'UniformOutput',false);
                    lambda_proj = cellfun(@(x,H) x.projection(H),lambda,basis,'UniformOutput',false);
                    lambda_ref_proj = cellfun(@(x,H) x.projection(H),lambda_ref,basis,'UniformOutput',false);
                    error_lambda_init = cellfun(@(x,xref) norm(x-xref)/norm(xref),lambda_proj,lambda_ref_proj,'UniformOutput',false);
                case 'residual'
                    error_U_init = 1;
                    error_w_init = repmat({1},[1,n]);
                    error_lambda_init = repmat({1},[1,n]);
            end
            
            % Check for convergence at initialization
            if error_U_init<=I.tolerance
                if I.display || I.displayIterations
                    fprintf('\nAlgorithm converged at initialization with error = %.3e w.r.t. U',error_U_init);
                    for k=1:n
                        fprintf('\n                                          error = %.3e w.r.t. w{%u}',error_w_init{k},k);
                        fprintf('\n                                          error = %.3e w.r.t. lambda{%u}',error_lambda_init{k},k);
                    end
                    fprintf('\n')
                    fprintf('\nElapsed time = %f s\n',toc(tTotal));
                end
                output.iteration = 0;
                output.errorGlobalSolutionInit = error_U_init;
                output.errorLocalSolutionInit = error_w_init;
                output.errorLagrangeMultiplierInit = error_lambda_init;
                return
            elseif I.displayIterations
                fprintf('\nInitialization : error = %.3e w.r.t. U',error_U_init);
                for k=1:n
                    fprintf('\n                 error = %.3e w.r.t. w{%u}',error_w_init{k},k);
                    fprintf('\n                 error = %.3e w.r.t. lambda{%u}',error_lambda_init{k},k);
                end
                fprintf('\n');
            end
            
            [error_U,error_w,error_lambda]...
                = initializeIndicators(I.maxIterations,patches);
            [stagnation_U,stagnation_w,stagnation_lambda]...
                = initializeIndicators(I.maxIterations,patches);
            time = initializeIndicators(I.maxIterations,patches);
            relaxationParameter = initializeIndicators(I.maxIterations,patches);
            % [sparsityRatio_U,sparsityRatio_w,sparsityRatio_lambda]...
            %     = initializeIndicators(I.maxIterations,patches);
            
            V = U;
            U_old = U;
            if ~isempty(glob.timeSolver)
                if glob.timeOrder>=1
                    vV = vU;
                end
                if glob.timeOrder>=2
                    aV = aU;
                end
            end
            
            % Relaxation parameter rho
            if ~I.updateRelaxationParameter
                rho = I.calcRelaxationParameterRandom(glob,patches,interfaces,s,bases,ls,rv,V,w);
            else
                rho = 1;
            end
            
            err_w = cell(1,n);
            err_lambda = cell(1,n);
            dim_w = cell(1,n);
            dim_lambda = cell(1,n);
            N = cell(1,n);
            y_w = cell(1,n);
            y_lambda = cell(1,n);
            
            % Iteration step - Main loop
            for iter=1:I.maxIterations
                
                tIter = tic;
                
                % Global problem
                % Global solution V
                indices = U.basis.indices;
                for k=1:n
                    indices = indices.addIndices(lambda{k}.basis.indices);
                end
                basis = SparseTensorProductFunctionalBasis(bases,indices);
                lambda_proj = cellfun(@(x) x.projection(basis),lambda,'UniformOutput',false);
                U_proj = U.projection(basis);
                V_old = V;
                V = V.projection(basis);
                b_out = glob.b_out;
                if isempty(glob.timeSolver)
                    % Time-independent problem
                    lambda_data = cellfun(@(x) x.data',lambda_proj,'UniformOutput',false);
                    U_data = U_proj.data';
                    V_data = V.data';
                    glob.b_out = b_out*[1,sparse(1,cardinal(basis)-1)];
                    V.data = glob.solve(interfaces,lambda_data,U_data,V_data)';
                    
%                     U_data = U_proj.data';
%                     V_data = V.data';
%                     for k=1:cardinal(basis)
%                         lambda_data = cellfun(@(x) x.data(k,:)',lambda_proj,'UniformOutput',false);
%                         if k~=1
%                             glob.b_out = [];
%                         end
%                         V.data(:,k) = glob.solve(interfaces,lambda_data,U_data(:,k),V_data(:,k));
%                     end
%                     V.data = V.data';
                else
                    % Time-dependent problem
                    T = gettimemodel(glob.timeSolver);
                    lambda_data = cell(cardinal(basis),1);
                    U_data = cell(cardinal(basis),1);
                    V_data = cell(cardinal(basis),1);
                    for k=1:cardinal(basis)
                        lambda_data{k} = cellfun(@(x) TIMEMATRIX(reshape(x.data(k,:,:),x.sz),T),lambda_proj,'UniformOutput',false);
                        U_data{k} = TIMEMATRIX(reshape(U_proj.data(k,:,:),U_proj.sz),T);
                        V_data{k} = TIMEMATRIX(reshape(V.data(k,:,:),V.sz),T);
                    end
                    if glob.timeOrder>=1
                        vU = vU.projection(basis);
                        vV = vV.projection(basis);
                        vU_data = cell(cardinal(basis),1);
                        vV_data = cell(cardinal(basis),1);
                        for k=1:cardinal(basis)
                            vU_data{k} = TIMEMATRIX(reshape(vU.data(k,:,:),vU.sz),T);
                            vV_data{k} = TIMEMATRIX(reshape(vV.data(k,:,:),vV.sz),T);
                        end
                    end
                    if glob.timeOrder>=2
                        aU = aU.projection(basis);
                        aV = aV.projection(basis);
                        aU_data = cell(cardinal(basis),1);
                        aV_data = cell(cardinal(basis),1);
                        for k=1:cardinal(basis)
                            aU_data{k} = TIMEMATRIX(reshape(aU.data(k,:,:),aU.sz),T);
                            aV_data{k} = TIMEMATRIX(reshape(aV.data(k,:,:),aV.sz),T);
                        end
                    end
                    V_datak = cell(cardinal(basis),1);
                    if glob.timeOrder>=1
                        vV_datak = cell(cardinal(basis),1);
                    end
                    if glob.timeOrder>=2
                        aV_datak = cell(cardinal(basis),1);
                    end
                    glob0 = glob;
                    glob0.b_out = [];
                    if glob.timeOrder==1
                        [V_datak{1},~,vV_datak{1}] = glob.solve(interfaces,lambda_data{1},U_data{1},V_data{1},vU_data{1},vV_data{1});
                        parfor k=2:cardinal(basis)
                            [V_datak{k},~,vV_datak{k}] = glob0.solve(interfaces,lambda_data{k},U_data{k},V_data{k},vU_data{k},vV_data{k});
                        end
                    elseif glob.timeOrder==2
                        [V_datak{1},~,vV_datak{1},aV_datak{1}] = glob.solve(interfaces,lambda_data{1},U_data{1},V_data{1},vU_data{1},vV_data{1},aU_data{1},aV_data{1});
                        parfor k=1:cardinal(basis)
                            [V_datak{k},~,vV_datak{k},aV_datak{k}] = glob0.solve(interfaces,lambda_data{k},U_data{k},V_data{k},vU_data{k},vV_data{k},aU_data{k},aV_data{k});
                        end
                    end
                    for k=1:cardinal(basis)
                        V.data(k,:,:) = V_datak{k};
                        if glob.timeOrder>=1
                            vV.data(k,:,:) = vV_datak{k};
                        end
                        if glob.timeOrder>=2
                            aV.data(k,:,:) = aV_datak{k};
                        end
                    end
                end
                glob.b_out = b_out;
                
                % Relaxation parameter rho
                if I.updateRelaxationParameter
                    rho = I.calcRelaxationParameterRandom(glob,patches,interfaces,s,bases,ls,rv,V,w,V_old,U,U_old,rho,iter);
                end
                relaxationParameter(iter) = rho;
                
                % Global iterate U
                U_old = U;
                U = U_proj; % U = U.projection(basis);
                U.data = rho*V.data + (1-rho)*U.data;
                if ~isempty(glob.timeSolver)
                    if glob.timeOrder>=1
                        vU.data = rho*vV.data + (1-rho)*vU.data;
                    end
                    if glob.timeOrder>=2
                        aU.data = rho*aV.data + (1-rho)*aU.data;
                    end
                end
                
                stagnation_U(iter) = norm(U-U_proj)/norm(U);
                
                switch lower(I.errorCriterion)
                    case 'none'
                        error_U(iter) = 1;
                    case 'reference'
                        if isempty(glob.timeSolver)
                            % Time-independent problem
                            U_out = U*glob.P_out';
                        else
                            % Time-dependent problem
                            U_out = U;
                            data = permute(U.data,[2 1 3]);
                            sz = size(data);
                            data = glob.P_out*data(:,:);
                            data = reshape(data,[size(data,1),sz(2:3)]);
                            U_out.data = permute(data,[2 1 3]);
                            sz_out = size(U_out.data);
                            U_out.sz = sz_out(2:end);
                        end
                        indices_ref = addIndices(U_out.basis.indices,U_ref.basis.indices);
                        basis_ref = SparseTensorProductFunctionalBasis(bases,indices_ref);
                        U_out_proj = U_out.projection(basis_ref);
                        U_ref_proj = U_ref.projection(basis_ref);
                        error_U(iter) = norm(U_out_proj-U_ref_proj)/norm(U_ref_proj);
                    case 'residual'
                        error_U(iter) = 1;
                end
                
                % sparsityRatio_U(iter) = getSparsityRatio(U.data);
                
                % Local problems
                % Local solutions (w{k},lambda{k}) without change of variable
                %                             z{k} with change of variable
                w_old = w;
                lambda_old = lambda;
                patch = patches.patches;
                interface = interfaces.interfaces;
                indices_init = U.basis.indices;
                parfor k=1:n
                    indices = indices_init;
                    indices = indices.addIndices(lambda{k}.basis.indices);
                    indices = indices.addIndices(w{k}.basis.indices);
                    basis = SparseTensorProductFunctionalBasis(bases,indices);
                    U_proj = U.projection(basis);
                    w_proj = w{k}.projection(basis);
                    lambda_proj = lambda{k}.projection(basis);
                    
                    if isempty(patch{k}.timeSolver)
                        fun = @(xi) solve(calcOperator(patch{k}.eval(xi)),interface{k},...
                            U_proj(xi)',w_proj(xi)',lambda_proj(xi)');
                    else
                        fun = @(xi) solve(calcOperator(patch{k}.eval(xi)),interface{k},...
                            TIMEMATRIX(reshape(U_proj(xi),U_proj.sz),gettimemodel(patch{k}.timeSolver)),...
                            TIMEMATRIX(reshape(w_proj(xi),w_proj.sz),gettimemodel(patch{k}.timeSolver)),...
                            TIMEMATRIX(reshape(lambda_proj(xi),lambda_proj.sz),gettimemodel(patch{k}.timeSolver)));
                    end
                    
                    fun = CellValuedUserDefinedFunction(fun,d,[2,1]);
                    fun.evaluationAtMultiplePoints = false;
                    
                    [u,err,x,y] = s.leastSquaresCell(fun,bases,ls,rv);
                    
                    [w{k},lambda{k}] = u{:};
                    [err_w{k}{iter},err_lambda{k}{iter}] = err{:};
                    [y_w{k}{iter},y_lambda{k}{iter}] = y{:};
                    dim_w{k}(iter) = cardinal(w{k}.basis);
                    dim_lambda{k}(iter) = cardinal(lambda{k}.basis);
                    N{k}(iter) = size(x,1);
                    
                    indices = addIndices(w{k}.basis.indices,w_old{k}.basis.indices);
                    basis = SparseTensorProductFunctionalBasis(bases,indices);
                    w_proj = w{k}.projection(basis);
                    w_old_proj = w_old{k}.projection(basis);
                    stagnation_w{k}(iter) = norm(w_proj-w_old_proj)/norm(w_proj);
                    
                    indices = addIndices(lambda{k}.basis.indices,lambda_old{k}.basis.indices);
                    basis = SparseTensorProductFunctionalBasis(bases,indices);
                    lambda_proj = lambda{k}.projection(basis);
                    lambda_old_proj = lambda_old{k}.projection(basis);
                    stagnation_lambda{k}(iter) = norm(lambda_proj-lambda_old_proj)/norm(lambda_proj);
                    
                    switch lower(I.errorCriterion)
                        case 'none'
                            error_w{k}(iter) = 1;
                            error_lambda{k}(iter) = 1;
                        case 'reference'
                            indices_ref = addIndices(w{k}.basis.indices,w_ref{k}.basis.indices);
                            basis_ref = SparseTensorProductFunctionalBasis(bases,indices_ref);
                            w_proj = w{k}.projection(basis_ref);
                            w_ref_proj = w_ref{k}.projection(basis_ref);
                            error_w{k}(iter) = norm(w_proj-w_ref_proj)/norm(w_ref_proj);
                            
                            indices_ref = addIndices(lambda{k}.basis.indices,lambda_ref{k}.basis.indices);
                            basis_ref = SparseTensorProductFunctionalBasis(bases,indices_ref);
                            lambda_proj = lambda{k}.projection(basis_ref);
                            lambda_ref_proj = lambda_ref{k}.projection(basis_ref);
                            error_lambda{k}(iter) = norm(lambda_proj-lambda_ref_proj)/norm(lambda_ref_proj);
                        case 'residual'
                            error_w{k}(iter) = 1;
                            error_w{k}(iter) = 1;
                    end
                    
                    % sparsityRatio_w{k}(iter) = getSparsityRatio(w{k}.data);
                    % sparsityRatio_lambda{k}(iter) = getSparsityRatio(lambda{k}.data);
                    
                end
                
                time(iter) = toc(tIter);
                
                % Check for convergence
                if error_U(iter)<=I.tolerance
                    break
                else
                    if I.displayIterations
                        fprintf('\nIteration #%2.d : stagnation = %.3e, error = %.3e w.r.t. U',iter,stagnation_U(iter),error_U(iter));
                        for k=1:n
                            fprintf('\n                stagnation = %.3e, error = %.3e w.r.t. w{%u}',stagnation_w{k}(iter),error_w{k}(iter),k);
                            fprintf('\n                stagnation = %.3e, error = %.3e w.r.t. lambda{%u}',stagnation_lambda{k}(iter),error_lambda{k}(iter),k);
                        end
                        fprintf('\n                elapsed time = %f s\n',time(iter));
                    end
                end
            end
            
            totalTime = toc(tTotal);
            
            if I.display || I.displayIterations
                if error_U(iter)<=I.tolerance
                    str = sprintf('Algorithm converged at iteration #%d with ',iter);
                    fprintf(['\n' str 'error = %.3e w.r.t. U'],error_U(iter));
                else
                    str = sprintf('Algorithm stopped at iteration #%d with ',iter);
                    fprintf(['\n' str 'error = %.3e w.r.t. U'],error_U(iter));
                end
                for k=1:n
                    fprintf(['\n' repmat(sprintf(' '),1,length(str)) 'error = %.3e w.r.t. w{%u}'],error_w{k}(iter),k);
                    fprintf(['\n' repmat(sprintf(' '),1,length(str)) 'error = %.3e w.r.t. lambda{%u}'],error_lambda{k}(iter),k);
                end
                fprintf('\n')
                fprintf('\nElapsed time = %f s\n',totalTime);
            end
            
            % Save outputs
            output.iteration = iter;
            output.totalTime = totalTime;
            output.time = time;
            output.relaxationParameter = relaxationParameter;
            
            output.errorGlobalSolutionInit = error_U_init;
            output.errorLocalSolutionInit = error_w_init;
            output.errorLagrangeMultiplierInit = error_lambda_init;
            
            output.errorGlobalSolution = error_U;
            output.errorLocalSolution = error_w;
            output.errorLagrangeMultiplier = error_lambda;
            
            output.stagnationGlobalSolution = stagnation_U;
            output.stagnationLocalSolution = stagnation_w;
            output.stagnationLagrangeMultiplier = stagnation_lambda;
            
            output.CVErrorLocalSolution = err_w;
            output.CVErrorLagrangeMultiplier = err_lambda;
            
            output.dimBasisLocalSolution = dim_w;
            output.dimBasisLagrangeMultiplier = dim_lambda;
            
            output.nbSamples = N;
            output.SamplesLocalSolution = y_w;
            output.SamplesLagrangeMultiplier = y_lambda;
            
        end

        function rho = calcRelaxationParameter(I,glob,patches,interfaces,V,w_old,V_old,U_old,U_old_old,rho_old,iter)
            % function rho = calcRelaxationParameter(I,glob,patches,interfaces,V,w_old,V_old,U_old,U_old_old,rho_old,iter)
            % Calculates relaxation parameter rho
            %
            % Inputs:
            % I: IterativeSolver
            % glob: Global
            % patches: Patches
            % interfaces: Interfaces
            % V: m_U-by-1 doubles containing current global solution V_{k}
            % w_old: 1-by-n cell of m_w-by-1 doubles containing previous local solution w_{k-1}
            % V_old: m_U-by-1 doubles containing previous global solution V_{k-1}
            % U_old: m_U-by-1 doubles containing previous global iterate U_{k-1}
            % U_old_old: m_U-by-1 doubles containing previous global iterate U_{k-2}
            % rho_old: 1-by-1 double containing previous relaxation parameter rho_{k-1}
            % iter: 1-by-1 double containing current iteration number k
            %
            % Outputs:
            % rho: 1-by-1 double containing current relaxation parameter rho_{k}
            
            if nargin<5 || (isempty(V) && isempty(w_old))
                [V,w_old,~] = initializeVariables(glob,patches,interfaces);
            elseif isempty(V)
                [V,~,~] = initializeVariables(glob,patches,interfaces);
            elseif nargin<6 || isempty(w_old)
                [~,w_old,~] = initializeVariables(glob,patches,interfaces);
            end
            if nargin<7
                [V_old,~,~] = initializeVariables(glob,patches,interfaces);
            end
            if nargin<8
                [U_old,~,~] = initializeVariables(glob,patches,interfaces);
            end
            if nargin<9
                [U_old_old,~,~] = initializeVariables(glob,patches,interfaces);
            end
            if nargin<10
                rho_old = 1;
            end
            if nargin<11
                iter = 0;
            end
            
            if isscalar(I.relaxation)
                rho = max(I.relaxation,eps);
                if I.display || I.displayIterations
                    fprintf('\nRelaxation parameter : rho = %f\n',rho);
                end
            elseif ischar(I.relaxation)
                switch lower(I.relaxation)
                    case {'optimal','optimalapproximation'}
                        t_rho = tic;
                        
                        A = @(W) calcIterationOperator(W,glob,patches,interfaces,V,w_old);
                        
                        if isempty(glob.timeSolver)
                            sz_W = getnbddlfree(glob.S);
                        else
                            sz_W = getnbddlfree(glob.S)*getnbtimedof(glob.timeSolver);
                        end
                        
                        eigopts = struct();
                        % eigopts.issym = false;
                        % eigopts.isreal = true;
                        eigopts.tol = eps;
                        if ~isempty(glob.solver) && isa(glob.solver,'NEWTONSOLVER')
                            tol = getparam(glob.solver,'tol');
                            eigopts.tol = max(tol,eigopts.tol);
                        end
                        n = numel(patches);
                        patch = patches.patches;
                        for k=1:n
                            if ~isempty(patch{k}.solver) && isa(patch{k}.solver,'NEWTONSOLVER')
                                tol = getparam(patch{k}.solver,'tol');
                                eigopts.tol = max(tol,eigopts.tol);
                            end
                        end
                        % eigopts.maxit = 300;
                        % eigopts.v0 = rand(sz_W,1);
                        % eigopts.disp = 0;
                        
                        [~,lambda_max,flag_max] = eigs(A,sz_W,1,'lm',eigopts);
                        
                        switch lower(I.relaxation)
                            case 'optimalapproximation'
                                rho = 1/lambda_max;
                                rho = max(rho,eps);
                                if I.display || I.displayIterations
                                    fprintf('\nApproximation of optimal relaxation parameter : rho = %f\n',rho);
                                    if flag_max
                                        fprintf('                                                lambda_max = %f, flag_max = %d -> eigenvalue not converged\n',lambda_max,flag_max);
                                    else
                                        fprintf('                                                lambda_max = %f, flag_max = %d -> eigenvalue converged\n',lambda_max,flag_max);
                                    end
                                    fprintf('                                                elapsed time = %f s\n',toc(t_rho));
                                end
                            case 'optimal'
                                [~,lambda_min,flag_min] = eigs(A,sz_W,1,'sm',eigopts);
                                rho = 2/(lambda_min + lambda_max);
                                rho = max(rho,eps);
                                if I.display || I.displayIterations
                                    fprintf('\nOptimal relaxation parameter : rho = %f\n',rho);
                                    if flag_max
                                        fprintf('                               lambda_max = %f, flag_max = %d -> eigenvalue not converged\n',lambda_max,flag_max);
                                    else
                                        fprintf('                               lambda_max = %f, flag_max = %d -> eigenvalue converged\n',lambda_max,flag_max);
                                    end
                                    if flag_min
                                        fprintf('                               lambda_min = %f, flag_min = %d -> eigenvalue not converged\n',lambda_min,flag_min);
                                    else
                                        fprintf('                               lambda_min = %f, flag_min = %d -> eigenvalue converged\n',lambda_min,flag_min);
                                    end
                                    fprintf('                               elapsed time = %f s\n',toc(t_rho));
                                end
                        end
                    case 'aitken'
                        if iter>=3
                            delta = V-U_old;
                            delta_old = V_old-U_old_old;
                            if isempty(glob.timeSolver)
                                rho = -rho_old*(delta-delta_old)'*delta_old/norm(delta-delta_old)^2;
                            else
                                rho = -rho_old*full(integratemtimes((delta-delta_old)',delta_old))/norm(delta-delta_old)^2;
                            end
                            rho = max(rho,eps);
                            if I.display || I.displayIterations
                                fprintf('\nRelaxation parameter dynamically updated through Aitken''s Delta-Squared method : rho = %f\n',rho);
                            end
                        else
                            rho = 1;
                            if I.display || I.displayIterations
                                fprintf('\nInitial relaxation parameter : rho = %f\n',rho);
                            end
                        end
                end
            else
                error('Relaxation should be a scalar or a char')
            end
            
        end
        
        function rho = calcRelaxationParameterRandom(I,glob,patches,interfaces,s,bases,ls,rv,V,w_old,V_old,U_old,U_old_old,rho_old,iter)
            % function rho = calcRelaxationParameterRandom(I,glob,patches,interfaces,s,bases,ls,rv,V,w_old,V_old,U_old,U_old_old,rho_old,iter)
            % Calculates relaxation parameter rho
            %
            % Inputs:
            % I: IterativeSolver
            % glob: Global
            % patches: Patches
            % interfaces: Interfaces
            % s: AdaptiveSparseTensorAlgorithm
            % bases: FunctionalBases
            % ls: LinearModelLearningSquareLoss
            % rv: RandomVector or RandomVariable (optional)
            % V: FunctionalBasisArray of size m_U of current global solution V_{k}
            % w_old: 1-by-n cell of FunctionalBasisArray of size m_w of previous local solution w_{k-1}
            % V_old: FunctionalBasisArray of size m_U of previous global solution V_{k-1}
            % U_old: FunctionalBasisArray of size m_U of previous global iterate U_{k-1}
            % U_old_old: FunctionalBasisArray of size m_U of previous global iterate U_{k-2}
            % rho_old: 1-by-1 double containing previous relaxation parameter rho_{k-1}
            % iter: 1-by-1 double containing current iteration number k
            %
            % Outputs:
            % rho: 1-by-1 double containing current relaxation parameter rho_{k}
            
            d = length(bases);
            rvb = getRandomVector(bases);
            if nargin<8 || isempty(rv)
                rv = rvb;
            elseif isa(rv,'RandomVariable')
                rv = RandomVector(rv,d);
            end
            if nargin<9 || (isempty(V) && isempty(w_old))
                [V,w_old,~] = initializeRandomVariables(glob,patches,interfaces,bases);
            elseif isempty(V)
                [V,~,~] = initializeRandomVariables(glob,patches,interfaces,bases);
            elseif nargin<10 || isempty(w_old)
                [~,w_old,~] = initializeRandomVariables(glob,patches,interfaces,bases);
            end
            if nargin<11
                [V_old,~,~] = initializeRandomVariables(glob,patches,interfaces,bases);
            end
            if nargin<12
                [U_old,~,~] = initializeRandomVariables(glob,patches,interfaces,bases);
            end
            if nargin<13
                [U_old_old,~,~] = initializeRandomVariables(glob,patches,interfaces,bases);
            end
            if nargin<14
                rho_old = 1;
            end
            if nargin<15
                iter = 0;
            end
            
            if isscalar(I.relaxation)
                rho = max(I.relaxation,eps);
                if I.display
                    fprintf('\nRelaxation parameter : rho = %f\n',rho);
                end
            elseif ischar(I.relaxation)
                switch lower(I.relaxation)
                    case {'optimal','optimalapproximation'}
                        t_rho = tic;
                        n = numel(patches);
                        
                        indices = V.basis.indices;
                        for k=1:n
                            indices = indices.addIndices(w_old{k}.basis.indices);
                        end
                        
                        A = @(W) calcIterationOperatorRandom(W,glob,patches,interfaces,s,bases,ls,rv,V,w_old);
                        
                        sz_W = getnbddlfree(glob.S)*cardinal(indices);
                        
                        eigopts = struct();
                        % eigopts.issym = false;
                        % eigopts.isreal = true;
                        eigopts.tol = s.tol;
                        if isfield(glob.solver,'tol')
                            eigopts.tol = max(glob.solver.tol,eigopts.tol);
                        end
                        patch = patches.patches;
                        for k=1:n
                            if isfield(patch{k}.solver,'tol')
                                eigopts.tol = max(patch{k}.solver.tol,eigopts.tol);
                            end
                        end
                        % eigopts.maxit = 300;
                        % eigopts.v0 = rand(sz_W,1);
                        % eigopts.disp = 0;
                        
                        [~,lambda_max,flag_max] = eigs(A,sz_W,1,'lm',eigopts);
                        % [~,lambda_max,flag_max] = power_iteration(A,sz_W,eigopts);
                        
                        switch lower(I.relaxation)
                            case 'optimalapproximation'
                                rho = 1/lambda_max;
                                rho = max(rho,eps);
                                if I.display
                                    fprintf('\nApproximation of optimal relaxation parameter : rho = %f\n',rho);
                                    if flag_max
                                        fprintf('                                                lambda_max = %f, flag_max = %d -> eigenvalue not converged\n',lambda_max,flag_max);
                                    else
                                        fprintf('                                                lambda_max = %f, flag_max = %d -> eigenvalue converged\n',lambda_max,flag_max);
                                    end
                                    fprintf('                                                elapsed time = %f s\n',toc(t_rho));
                                end
                            case 'optimal'
                                [~,lambda_min,flag_min] = eigs(A,sz_W,1,'sm',eigopts);
                                rho = 2/(lambda_min + lambda_max);
                                rho = max(rho,eps);
                                if I.display
                                    fprintf('\nOptimal relaxation parameter : rho = %f\n',rho);
                                    if flag_max
                                        fprintf('                               lambda_max = %f, flag_max = %d -> eigenvalue not converged\n',lambda_max,flag_max);
                                    else
                                        fprintf('                               lambda_max = %f, flag_max = %d -> eigenvalue converged\n',lambda_max,flag_max);
                                    end
                                    if flag_min
                                        fprintf('                               lambda_min = %f, flag_min = %d -> eigenvalue not converged\n',lambda_min,flag_min);
                                    else
                                        fprintf('                               lambda_min = %f, flag_min = %d -> eigenvalue converged\n',lambda_min,flag_min);
                                    end
                                    fprintf('                               elapsed time = %f s\n',toc(t_rho));
                                end
                        end
                    case 'aitken'
                        if iter>=3
                            indices = addIndices(V.basis.indices,addIndices(U_old.basis.indices,U_old_old.basis.indices));
                            basis = SparseTensorProductFunctionalBasis(bases,indices);
                            V_proj = V.projection(basis);
                            U_old_proj = U_old.projection(basis);
                            delta = V_proj-U_old_proj;
                            V_old_proj = V_old.projection(basis);
                            U_old_old_proj = U_old_old.projection(basis);
                            delta_old = V_old_proj-U_old_old_proj;
                            rho = -rho_old*full(sum(sum(dot(delta-delta_old,delta_old))))/norm(delta-delta_old)^2;
                            % rho = -rho_old*full(sum(sum((delta.data(:,:)-delta_old.data(:,:)).*delta_old.data(:,:))))/norm(delta-delta_old)^2;
                            rho = max(rho,eps);
                            fprintf('\nRelaxation parameter dynamically updated through Aitken''s Delta-Squared method : rho = %f\n',rho);
                        else
                            rho = 1;
                            fprintf('\nInitial relaxation parameter : rho = %f\n',rho);
                        end
                end
            else
                error('Relaxation parameter rho should be a scalar or a char')
            end
            
        end
        
        function [glob,patches] = initializeRightHandSide(I,glob,patches)
            % function [glob,patches] = initializeRightHandSide(I,glob,patches)
            
            n = numel(patches);
            
            if ~isempty(I.timeSolver)
                T = gettimemodel(I.timeSolver);
            end
            if isempty(glob.b)
                sz_U = getnbddlfree(glob.S);
                glob.b = sparse(sz_U,1);
                if ~isempty(I.timeSolver)
                    glob.b = glob.b*zero(T);
                end
            end
            for k=1:n
                if isempty(patches.patches{k}.b)
                    sz_w = getnbddl(patches.patches{k}.S);
                    patches.patches{k}.b = sparse(sz_w,1);
                    if ~isempty(I.timeSolver)
                        patches.patches{k}.b = patches.patches{k}.b*zero(T);
                    end
                end
            end
        end

    end
    
end
