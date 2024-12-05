classdef DirectSolver
    %
    
    %
    %   A MATLAB class for representing a multiscale direct solver
    %
    
    properties
        changeOfVariable
        solver
        initializationType
        timeSolver
        timeOrder
        display
    end
    
    methods
        
        % Constructor
        function D = DirectSolver(varargin)
            % Class DirectSolver
            %
            % function D = DirectSolver(varargin)
            % D.changeOfVariable: reformulation of local problems as single-field problems by introducing a change of variable w=w_tilde+z (true or false), false by default
            % D.solver: solver for multiscale problem ([] or NEWTONSOLVER), [] by default
            % if D.solver is an instance of class NEWTONSOLVER,
            % D.initializationType: type of initialization for iterative resolution of multiscale problem ('zero' or 'one'), 'zero' by default
            % D.timeSolver: time solver for multiscale problem ([], DGTIMESOLVER, EULERSOLVER, NEWMARKSOLVER), [] by default
            % D.timeOrder: time order for multiscale problem, 0 by default
            % D.display: display (true or false), true by default
            
            expectedInitializationTypes = {'zero','one'};
            
            p = ImprovedInputParser;
            addParameter(p,'changeOfVariable',false,@islogical);
            addParameter(p,'solver',[]);
            % addParameter(p,'initializationType','zero',@ischar);
            addParameter(p,'initializationType','zero',...
                @(x) any(validatestring(x,expectedInitializationTypes)))
            addParameter(p,'timeSolver',[]);
            addParameter(p,'timeOrder',0,@isscalar);
            addParameter(p,'display',true,@islogical);
            parse(p,varargin{:});
            D = passMatchedArgsToProperties(p,D);
            
        end
        
        function [U,w,lambda,output,vU,vw,vlambda,aU,aw,alambda] = solve(D,glob,patches,interfaces,varargin)
            % function [U,w,lambda,output] = solve(D,glob,patches,interfaces)
            % function [U,w,lambda,output,vU,vw,vlambda] = solve(D,glob,patches,interfaces)
            % function [U,w,lambda,output,vU,vw,vlambda,aU,aw,alambda] = solve(D,glob,patches,interfaces)
            % function u = solve(D,glob,patches,interfaces)
            % Solves multiscale deterministic problem based on
            % non-overlapping domain decomposition using direct solver
            % If the number of output arguments is equal to one, it returns
            % 1-by-(1+2xn) cell representing multiscale solution u and containing 
            % both global solution U, local solution w and  Lagrange multiplier lambda
            %
            % Inputs:
            % D: DirectSolver
            % glob: GlobalOutside
            % patches: Patches
            % interfaces: Interfaces
            %
            % Outputs:
            % U: m_U-by-p doubles containing global solution U
            % w: 1-by-n cell of m_w-by-p doubles containing local solution w
            % lambda: 1-by-n cell of m_l-by-p doubles containing Lagrange multiplier lambda
            % u: (1+2xn)-by-1 cell of doubles arrays containing multiscale solution u=(U,w,lambda)
            % output.time: 1-by-1 double containing CPU time t
            % For time-dependent problems,
            % output.result: structure containing outputs of time solver
            % For first- and second-order time-dependent problems,
            % vU: m_U-by-p doubles containing global velocity vU
            % vw: 1-by-n cell of m_w-by-p doubles containing local velocity vw
            % vlambda: 1-by-n cell of m_l-by-p doubles containing Lagrange multiplier velocity vlambda
            % v: (1+2xn)-by-1 cell of doubles arrays containing multiscale velocity v=(vU,vw,vlambda)
            % For second-order time-dependent problems,
            % aU: m_U-by-p doubles containing global acceleration vU
            % aw: 1-by-n cell of m_w-by-p doubles containing local acceleration aw
            % alambda: 1-by-n cell of m_l-by-p doubles containing Lagrange multiplier acceleration alambda
            % a: (1+2xn)-by-1 cell of doubles arrays containing multiscale acceleration a=(aU,aw,alambda)
            % where
            % n is the number of patches
            % m_U is the dimension of the spatial approximation space of global solution U
            % m_w is the dimension of the spatial approximation space of local solution w
            % m_l is the dimension of the spatial approximation space of Lagrange multiplier lambda
            % p is the dimension of the time approximation space of multiscale problem
            
            t = tic;
            
            if D.display
                fprintf('\n ---------------------------------------');
                fprintf('\n ------------ Direct solver ------------')
                fprintf('\n ---------------------------------------\n');
            end
            
            [glob,patches] = initializeRightHandSide(D,glob,patches);
            
            if D.changeOfVariable
                patches = setBoundaryCondition(patches,interfaces);
            end  
                
            % Stiffness matrix A
            [A,Atang] = assembleStiffnessOperator(D,glob,patches,interfaces);
            
            if ~isempty(D.timeSolver)
                % Mass matrix M
                M = assembleMassOperator(D,glob,patches,interfaces);
                if D.timeOrder==2
                    % Damping matrix C
                    C = assembleDampingOperator(D,glob,patches,interfaces);
                end
            end
            
            % Sollicitation vector b
            b = assembleRightHandSide(D,glob,patches,interfaces);
            
            % Multiscale solution u
            if isempty(D.timeSolver)
                % Time-independent problem
                if isempty(D.solver)
                    u = A\b;
                else
                    if D.changeOfVariable && isanlsolver(D.solver)
                        warning('Wrong local transfert from patch to interface. Set parameter ''changeOfVariable'' to false.')
                    end
                    if isa(D.solver,'NEWTONSOLVER')
                        u0 = setInitialSolution(D.initializationType,b);
                        u = solve(D.solver,b,A,Atang,u0);
                    end
                end
            else
                % Time-dependent problem
                if D.timeOrder==1
                    u0 = setInitialCondition(D,glob,patches,interfaces);
                    [u,output.result,v] = dsolve(D.timeSolver,b,M,A,u0);
                elseif D.timeOrder==2
                    [u0,v0] = setInitialCondition(D,glob,patches,interfaces);
                    [u,output.result,v,a] = ddsolve(D.timeSolver,b,M,A,C,u0,v0);
                end
            end
            
            if ~D.changeOfVariable
                % Without change of variable
                % Multiscale solution u=(U,w{k},lambda{k})
                [u,U,w,lambda] = getSolution(D,glob,patches,interfaces,u);
                if ~isempty(D.timeSolver)
                    if D.timeOrder>=1
                        [v,vU,vw,vlambda] = getSolution(D,glob,patches,interfaces,v);
                    end
                    if D.timeOrder>=2
                        [a,aU,aw,alambda] = getSolution(D,glob,patches,interfaces,a);
                    end
                end
            else
                % With change of variable
                % Multiscale solution u=(U,z{k})
                % Local solution w{k}
                % Lagrange multiplier lambda{k}
                [u,U,w] = getSolution(D,glob,patches,interfaces,u);
                if isempty(D.timeSolver)
                    lambda = computeLagrangeMultiplier(patches,interfaces,w);
                else
                    [v,vU,vw] = getSolution(D,glob,patches,interfaces,v);
                    if D.timeOrder==1
                        lambda = computeLagrangeMultiplier(patches,interfaces,w,vw);
                    elseif D.timeOrder==2
                        [a,aU,aw] = getSolution(D,glob,patches,interfaces,a);
                        lambda = computeLagrangeMultiplier(patches,interfaces,w,vw,aw);
                    end
                    T = gettimemodel(D.timeSolver);
                    tSolver = EULERTIMESOLVER(T,'eulertype','implicit','display',false);
                    if D.timeOrder>=1
                        vlambda = cellfun(@(lambda) derivative(tSolver,lambda),lambda,'UniformOutput',false);
                    end
                    if D.timeOrder>=2
                        alambda = cellfun(@(vlambda) derivative(tSolver,vlambda),vlambda,'UniformOutput',false);
                    end
                end
                    
                % Multiscale solution u=(U,w{k},lambda{k})
                u = [w;lambda];
                u = [{U};u(:)];
                if ~isempty(D.timeSolver)
                    u = cellfun(@(u) getvalue(u),u,'UniformOutput',false);
                    if D.timeOrder>=1
                        % Multiscale velocity v=(vU,vw{k},vlambda{k})
                        v = [vw;vlambda];
                        v = [{vU};v(:)];
                        v = cellfun(@(v) getvalue(v),v,'UniformOutput',false);
                    end
                    if D.timeOrder>=2
                        % Multiscale acceleration a=(aU,aw{k},alambda{k})
                        a = [aw;alambda];
                        a = [{aU};a(:)];
                        a = cellfun(@(a) getvalue(a),a,'UniformOutput',false);
                    end
                end
            end
            
            output.time = toc(t);
            
            if D.display
                fprintf('\nElapsed time = %f s\n',output.time);
            end
            
            if nargout==1
                % U = vertcat(u{:});
                U = u;
            end
            
        end
        
        function [U,w,lambda,output] = solveRandom(D,glob,patches,interfaces,s,bases,ls,rv)
            % function [U,w,lambda,output] = solveRandom(D,glob,patches,interfaces,s,bases,ls,rv)
            % Solves multiscale stochastic problem based on non-overlapping
            % domain decomposition using direct solver
            %
            % Inputs:
            % D: DirectSolver
            % glob: GlobalOutside
            % patches: Patches
            % interfaces: Interfaces
            % s: AdaptiveSparseTensorAlgorithm
            % bases: FunctionalBases
            % ls: LinearModelLearningSquareLoss
            % rv: RandomVector or RandomVariable (optional)
            %
            % Outputs:
            % U: FunctionalBasisArray of size m_U of global solution U
            % w: 1-by-n cell of FunctionalBasisArray of size m_w of local solution w
            % lambda: 1-by-n cell of FunctionalBasisArray of size m_l of Lagrange multiplier lambda
            % output.time: 1-by-1 double containing CPU time t
            % output.CVErrorGlobalSolution: 1-by-m_U doubles containing the cross-validation error of global solution U
            % output.CVErrorLocalSolution: 1-by-n cell of 1-by-m_w doubles containing the cross-validation error of local solution w
            % output.CVErrorLagrangeMultiplier: 1-by-n cell of 1-by-m_l doubles containing the cross-validation error of Lagrange mutliplier lambda
            % output.nbSamples: 1-by-1 double containing sample size N
            % output.SamplesGlobalSolution: N-by-m_U doubles containing the evaluations of global solution U
            % output.SamplesLocalSolution: 1-by-n cell of N-by-m_w doubles containing the evaluations of local solution w
            % output.SamplesLagrangeMultiplier: 1-by-n cell of N-by-m_l doubles containing the evaluations of Lagrange multiplier lambda
            % where
            % n is the number of patches
            % N is the number of samples
            % m_U is the dimension of the spatial approximation space of global solution U
            % m_w is the dimension of the spatial approximation space of local solution w
            % m_l is the dimension of the spatial approximation space of Lagrange multiplier lambda
            
            d = length(bases);
            rvb = getRandomVector(bases);
            if nargin<8 || isempty(rv)
                rv = rvb;
            elseif isa(rv,'RandomVariable')
                rv = RandomVector(rv,d);
            end
            
            n = numel(patches);
            
            t = tic;
            
            if D.display
                fprintf('\n ---------------------------------------');
                fprintf('\n ------------ Direct solver ------------')
                fprintf('\n ---------------------------------------\n');
            end
            
            display_ = D.display;
            D.display = false;
            
            fun = @(xi) D.solve(glob,calcOperator(patches.eval(xi)),interfaces);
            fun = CellValuedUserDefinedFunction(fun,d,[1+2*n,1]);
            fun.evaluationAtMultiplePoints = false;
            
            [u,err,x,y] = s.leastSquaresCell(fun,bases,ls,rv);
            
            patch = patches.patches;
            interface = interfaces.interfaces;
            
            U = u{1};
            w = cellfun(@(patch) u{2*patch.number},patch,'UniformOutput',false);
            lambda = cellfun(@(interface) u{2*interface.number+1},interface,'UniformOutput',false);
            output.CVErrorGlobalSolution = err{1};
            output.CVErrorLocalSolution = cellfun(@(patch) err{2*patch.number},patch,'UniformOutput',false);
            output.CVErrorLagrangeMultiplier = cellfun(@(interface) err{2*interface.number+1},interface,'UniformOutput',false);
            output.nbSamples = size(x,1);
            output.SamplesGlobalSolution = y{1};
            output.SamplesLocalSolution = cellfun(@(patch) y{2*patch.number},patch,'UniformOutput',false);
            output.SamplesLagrangeMultiplier = cellfun(@(interface) y{2*interface.number+1},interface,'UniformOutput',false);
            
            output.time = toc(t);
            
            if display_
                fprintf('\nElapsed time = %f s\n',output.time);
            end
            
        end
        
        function [glob,patches] = initializeRightHandSide(D,glob,patches)
            % function [glob,patches] = initializeRightHandSide(D,glob,patches)
            
            n = numel(patches);
            
            if ~isempty(D.timeSolver)
                T = gettimemodel(D.timeSolver);
            end
            if isempty(glob.b)
                sz_U = getnbddlfree(glob.S);
                glob.b = sparse(sz_U,1);
                if ~isempty(D.timeSolver)
                    glob.b = glob.b*zero(T);
                end
            end
            for k=1:n
                if isempty(patches.patches{k}.b)
                    sz_w = getnbddl(patches.patches{k}.S);
                    patches.patches{k}.b = sparse(sz_w,1);
                    if ~isempty(D.timeSolver)
                        patches.patches{k}.b = patches.patches{k}.b*zero(T);
                    end
                end
            end
        end
        
        function b = assembleRightHandSide(D,glob,patches,interfaces)
            % function b = assembleRightHandSide(D,glob,patches,interfaces)
            
            n = numel(patches);
            patch = patches.patches;
            interface = interfaces.interfaces;
            
            if ~D.changeOfVariable
                % Without change of variable
                sz_lambda = cellfun(@(interface) getnbddl(interface.S),interface,'UniformOutput',false);
                
                b = cell(1+2*n,1);
                b{1} = glob.b;
                for k=1:n
                    b{2*k} = patch{k}.b;
                    b{2*k+1} = sparse(sz_lambda{k},1);
                    if ~isempty(D.timeSolver)
                        T = gettimemodel(D.timeSolver);
                        b{2*k+1} = b{2*k+1}*zero(T);
                    end
                end
            else
                % With change of variable
                P = cellfun(@(interface) interface.P_patch'*interface.P_globOut,interface,'UniformOutput',false);
                
                b = cell(1+n,1);
                b{1} = glob.b;
                for k=1:n
                    b{1} = b{1} + P{k}'*patch{k}.b;
                    b{k+1} = freevector(patch{k}.S,patch{k}.b);
                end
            end
            b = vertcat(b{:});
        end
        
        function [A,Atang] = assembleStiffnessOperator(D,glob,patches,interfaces)
            % function [A,Atang] = assembleStiffnessOperator(D,glob,patches,interfaces)
            
            n = numel(patches);
            patch = patches.patches;
            interface = interfaces.interfaces;
            
            if isempty(D.solver)
                if ~D.changeOfVariable
                    % Without change of variable
                    sz_U = getnbddlfree(glob.S);
                    sz_w = cellfun(@(patch) getnbddl(patch.S),patch,'UniformOutput',false);
                    sz_lambda = cellfun(@(interface) getnbddl(interface.S),interface,'UniformOutput',false);
                    
                    B_globOut = cellfun(@(interface) interface.P_globOut'*interface.M,interface,'UniformOutput',false);
                    B_patch = cellfun(@(interface) interface.P_patch'*interface.M,interface,'UniformOutput',false);
                    
                    A = cell(1+2*n,1+2*n);
                    A{1,1} = glob.A;
                    for k=1:n
                        for j=1:n
                            A{2*k,2*j} = sparse(sz_w{k},sz_w{j});
                            A{2*k+1,2*j+1} = sparse(sz_lambda{k},sz_lambda{j});
                            A{2*k,2*j+1} = sparse(sz_w{k},sz_lambda{j});
                            A{2*k+1,2*j} = sparse(sz_lambda{k},sz_w{j});
                        end
                        A{2*k,2*k} = patch{k}.A;
                        A{1,2*k} = sparse(sz_U,sz_w{k});
                        A{1,2*k+1} = B_globOut{k};
                        A{2*k,2*k+1} = -B_patch{k};
                        A{2*k,1} = A{1,2*k}';
                        A{2*k+1,1} = A{1,2*k+1}';
                        A{2*k+1,2*k} = A{2*k,2*k+1}';
                    end
                else
                    % With change of variable
                    sz_z = cellfun(@(patch) getnbddlfree(patch.S),patch,'UniformOutput',false);
                    
                    P = cellfun(@(interface) interface.P_patch'*interface.P_globOut,interface,'UniformOutput',false);
                    
                    A = cell(1+n,1+n);
                    A{1,1} = glob.A;
                    for k=1:n
                        for j=1:n
                            A{k+1,j+1} = sparse(sz_z{k},sz_z{j});
                        end
                        A{1,1} = A{1,1} + P{k}'*patch{k}.A*P{k};
                        A{k+1,k+1} = freematrix(patch{k}.S,patch{k}.A);
                        A{1,k+1} = P{k}'*freevector(patch{k}.S,patch{k}.A,2);
                        A{k+1,1} = freevector(patch{k}.S,patch{k}.A,1)*P{k};
                    end
                end
                A = cell2mat(A);
                Atang = [];
            else
                if isa(D.solver,'NEWTONSOLVER')
                    A = @(u) calcStiffnessOperator(D,glob,patches,interfaces,u);
                    Atang = @(u) calcTangentStiffnessOperator(D,glob,patches,interfaces,u);
                end
            end
            
        end
        
        function A = calcStiffnessOperator(D,glob,patches,interfaces,u)
            % function A = calcStiffnessOperator(D,glob,patches,interfaces,u)
            
            n = numel(patches);
            
            patch = patches.patches;
            interface = interfaces.interfaces;
            
            sz_U = getnbddlfree(glob.S);
            if ~D.changeOfVariable
                % Without change of variable
                sz_w = cellfun(@(patch) getnbddl(patch.S),patch,'UniformOutput',false);
                sz_lambda = cellfun(@(interface) getnbddl(interface.S),interface,'UniformOutput',false);
                rep = [sz_w{:};sz_lambda{:}];
                rep = [sz_U;rep(:)];
                u = mat2cell(u,rep);
                U = u{1};
                w = cellfun(@(patch) u{2*patch.number},patch,'UniformOutput',false);
                lambda = cellfun(@(interface) u{2*interface.number+1},interface,'UniformOutput',false);
                
                A = cell(1+2*n,1);
                if isa(glob.A,'function_handle')
                    A{1} = glob.A(U);
                else
                    A{1} = glob.A*U;
                end
                B_globOut = cellfun(@(interface) interface.P_globOut'*interface.M,interface,'UniformOutput',false);
                B_patch = cellfun(@(interface) interface.P_patch'*interface.M,interface,'UniformOutput',false);
                for k=1:n
                    A{1} = A{1} + B_globOut{k}*lambda{k};
                    if isa(patch{k}.A,'function_handle')
                        A{2*k} = patch{k}.A(w{k}) - B_patch{k}*lambda{k};
                    else
                        A{2*k} = patch{k}.A*w{k} - B_patch{k}*lambda{k};
                    end
                    A{2*k+1} = B_globOut{k}'*U - B_patch{k}'*w{k};
                end
            else
                % With change of variable
                sz_z = cellfun(@(patch) getnbddlfree(patch.S),patch,'UniformOutput',false);
                rep = [sz_z{:}];
                rep = [sz_U;rep(:)];
                u = mat2cell(u,rep);
                U = u{1};
                z = cellfun(@(patch) u{1+patch.number},patch,'UniformOutput',false);
                
                A = cell(1+n,1);
                if isa(glob.A,'function_handle')
                    A{1} = glob.A(U);
                else
                    A{1} = glob.A*U;
                end
                P = cellfun(@(interface) interface.P_patch'*interface.P_globOut,interface,'UniformOutput',false);
                for k=1:n
                    if isa(patch{k}.A,'function_handle')
                        A{1} = A{1} + P{k}'*patch{k}.A(P{k}*U) + P{k}'*patch{k}.A(z{k});
                        A{k+1} = freevector(patch{k}.S,patch{k}.A(P{k}*U),1) + freevector(patch{k}.S,patch{k}.A(z{k}));
                    else
                        A{1} = A{1} + P{k}'*patch{k}.A*P{k}*U + P{k}'*freevector(patch{k}.S,patch{k}.A,2)*z{k};
                        A{k+1} = freevector(patch{k}.S,patch{k}.A,1)*P{k}*U + freematrix(patch{k}.S,patch{k}.A)*z{k};
                    end
                end
            end
            A = cell2mat(A);
        end
        
        function Atang = calcTangentStiffnessOperator(D,glob,patches,interfaces,u)
            % function Atang = calcTangentStiffnessOperator(D,glob,patches,interfaces,u)
            
            n = numel(patches);
            
            patch = patches.patches;
            interface = interfaces.interfaces;
            
            sz_U = getnbddlfree(glob.S);
            if ~D.changeOfVariable
                % Without change of variable
                sz_w = cellfun(@(patch) getnbddl(patch.S),patch,'UniformOutput',false);
                sz_lambda = cellfun(@(interface) getnbddl(interface.S),interface,'UniformOutput',false);
                rep = [sz_w{:};sz_lambda{:}];
                rep = [sz_U;rep(:)];
                u = mat2cell(u,rep);
                U = u{1};
                w = cellfun(@(patch) u{2*patch.number},patch,'UniformOutput',false);
                
                Atang = cell(1+2*n,1+2*n);
                if ~isempty(glob.Atang) && isa(glob.Atang,'function_handle')
                    Atang{1,1} = glob.Atang(U);
                else
                    Atang{1,1} = glob.A;
                end
                for k=1:n
                    for j=1:n
                        Atang{2*k,2*j} = sparse(sz_w{k},sz_w{j});
                        Atang{2*k+1,2*j+1} = sparse(sz_lambda{k},sz_lambda{j});
                        Atang{2*k,2*j+1} = sparse(sz_w{k},sz_lambda{j});
                        Atang{2*k+1,2*j} = sparse(sz_lambda{k},sz_w{j});
                    end
                    B_globOut = interface{k}.P_globOut'*interface{k}.M;
                    B_patch = interface{k}.P_patch'*interface{k}.M;
                    if ~isempty(patch{k}.Atang) && isa(patch{k}.Atang,'function_handle')
                        Atang{2*k,2*k} = patch{k}.Atang(w{k});
                    else
                        Atang{2*k,2*k} = patch{k}.A;
                    end
                    Atang{1,2*k} = sparse(sz_U,sz_w{k});
                    Atang{1,2*k+1} = B_globOut;
                    Atang{2*k,2*k+1} = -B_patch;
                    Atang{2*k,1} = Atang{1,2*k}';
                    Atang{2*k+1,1} = Atang{1,2*k+1}';
                    Atang{2*k+1,2*k} = Atang{2*k,2*k+1}';
                end
            else
                % With change of variable
                sz_z = cellfun(@(patch) getnbddlfree(patch.S),patch,'UniformOutput',false);
                rep = [sz_z{:}];
                rep = [sz_U;rep(:)];
                u = mat2cell(u,rep);
                U = u{1};
                z = cellfun(@(patch) u{1+patch.number},patch,'UniformOutput',false);
                    
                Atang = cell(1+n,1+n);
                if ~isempty(glob.Atang) && isa(glob.Atang,'function_handle')
                    Atang{1,1} = glob.Atang(U);
                else
                    Atang{1,1} = glob.A;
                end
                for k=1:n
                    for j=1:n
                        Atang{k+1,j+1} = sparse(sz_z{k},sz_z{j});
                    end
                    P = interface{k}.P_patch'*interface{k}.P_globOut;
                    if ~isempty(patch{k}.Atang) && isa(patch{k}.Atang,'function_handle')
                        Atang{1,1} = Atang{1,1} + P'*patch{k}.Atang(P*U)*P;
                        Atang{k+1,k+1} = freematrix(patch{k}.S,patch{k}.Atang(z{k}));
                        Atang{1,k+1} = P'*freevector(patch{k}.S,patch{k}.Atang(z{k}),2);
                        Atang{k+1,1} = freevector(patch{k}.S,patch{k}.Atang(P*U),1)*P;
                    else
                        Atang{1,1} = Atang{1,1} + P'*patch{k}.A*P;
                        Atang{k+1,k+1} = freematrix(patch{k}.S,patch{k}.A);
                        Atang{1,k+1} = P'*freevector(patch{k}.S,patch{k}.A,2);
                        Atang{k+1,1} = freevector(patch{k}.S,patch{k}.A,1)*P;
                    end
                end
            end
            Atang = cell2mat(Atang);
        end
        
        function M = assembleMassOperator(D,glob,patches,interfaces)
            % function M = assembleMassOperator(D,glob,patches,interfaces)
            
            n = numel(patches);
            patch = patches.patches;
            interface = interfaces.interfaces;
            
            if ~D.changeOfVariable
                % Without change of variable
                sz_U = getnbddlfree(glob.S);
                sz_w = cellfun(@(patch) getnbddl(patch.S),patch,'UniformOutput',false);
                sz_lambda = cellfun(@(interface) getnbddl(interface.S),interface,'UniformOutput',false);
                
                M = cell(1+2*n,1+2*n);
                M{1,1} = glob.M;
                for k=1:n
                    for j=1:n
                        M{2*k,2*j} = sparse(sz_w{k},sz_w{j});
                        M{2*k+1,2*j+1} = sparse(sz_lambda{k},sz_lambda{j});
                        M{2*k,2*j+1} = sparse(sz_w{k},sz_lambda{j});
                        M{2*k+1,2*j} = sparse(sz_lambda{k},sz_w{j});
                    end
                    M{2*k,2*k} = patch{k}.M;
                    M{1,2*k} = sparse(sz_U,sz_w{k});
                    M{1,2*k+1} = sparse(sz_U,sz_lambda{k});
                    M{2*k,1} = M{1,2*k}';
                    M{2*k+1,1} = M{1,2*k+1}';
                end
            else
                % With change of variable
                sz_z = cellfun(@(patch) getnbddlfree(patch.S),patch,'UniformOutput',false);
                
                P = cellfun(@(interface) interface.P_patch'*interface.P_globOut,interface,'UniformOutput',false);
                
                M = cell(1+n,1+n);
                M{1,1} = glob.M;
                for k=1:n
                    for j=1:n
                        M{k+1,j+1} = sparse(sz_z{k},sz_z{j});
                    end
                    M{1,1} = M{1,1} + P{k}'*patch{k}.M*P{k};
                    M{k+1,k+1} = freematrix(patch{k}.S,patch{k}.M);
                    M{1,k+1} = P{k}'*freevector(patch{k}.S,patch{k}.M,2);
                    M{k+1,1} = freevector(patch{k}.S,patch{k}.M,1)*P{k};
                end
            end
            M = cell2mat(M);
        end
        
        function C = assembleDampingOperator(D,glob,patches,interfaces)
            % function C = assembleDampingOperator(D,glob,patches,interfaces)
            
            n = numel(patches);
            patch = patches.patches;
            interface = interfaces.interfaces;
            
            sz_U = getnbddlfree(glob.S);
            if ~D.changeOfVariable
                % Without change of variable
                sz_w = cellfun(@(patch) getnbddl(patch.S),patch,'UniformOutput',false);
                sz_lambda = cellfun(@(interface) getnbddl(interface.S),interface,'UniformOutput',false);
                
                C = cell(1+2*n,1+2*n);
                if ~isempty(glob.C)
                    C{1,1} = glob.C;
                else
                    C{1,1} = sparse(sz_U,sz_U);
                end
                for k=1:n
                    for j=1:n
                        C{2*k,2*j} = sparse(sz_w{k},sz_w{j});
                        C{2*k+1,2*j+1} = sparse(sz_lambda{k},sz_lambda{j});
                        C{2*k,2*j+1} = sparse(sz_w{k},sz_lambda{j});
                        C{2*k+1,2*j} = sparse(sz_lambda{k},sz_w{j});
                    end
                    if ~isempty(patch{k}.C)
                        C{2*k,2*k} = patch{k}.C;
                    else
                        C{2*k,2*k} = sparse(sz_w{k},sz_w{k});
                    end
                    C{1,2*k} = sparse(sz_U,sz_w{k});
                    C{1,2*k+1} = sparse(sz_U,sz_lambda{k});
                    C{2*k,1} = C{1,2*k}';
                    C{2*k+1,1} = C{1,2*k+1}';
                end
            else
                % With change of variable
                sz_z = cellfun(@(patch) getnbddlfree(patch.S),patch,'UniformOutput',false);
                
                P = cellfun(@(interface) interface.P_patch'*interface.P_globOut,interface,'UniformOutput',false);
                
                C = cell(1+n,1+n);
                if ~isempty(glob.C)
                    C{1,1} = glob.C;
                else
                    C{1,1} = sparse(sz_U,sz_U);
                end
                for k=1:n
                    for j=1:n
                        C{k+1,j+1} = sparse(sz_z{k},sz_z{j});
                    end
                    if ~isempty(patch{k}.C)
                        C{1,1} = C{1,1} + P{k}'*patch{k}.C*P{k};
                        C{k+1,k+1} = freematrix(patch{k}.S,patch{k}.C);
                        C{1,k+1} = P{k}'*freevector(patch{k}.S,patch{k}.C,2);
                        C{k+1,1} = freevector(patch{k}.S,patch{k}.C,1)*P{k};
                    else
                        C{1,k+1} = sparse(sz_U,sz_z{k});
                        C{k+1,1} = C{1,k+1}';
                    end
                end
            end
            C = cell2mat(C);
        end
        
        function [u0,v0] = setInitialCondition(D,glob,patches,interfaces)
            % function u0 = setInitialCondition(D,glob,patches,interfaces)
            % Sets multiscale initial solution vector u0
            % function [u0,v0] = setInitialCondition(D,glob,patches,interfaces)
            % Sets multiscale initial solution vector u0 and velocity vector v0
            %
            % Inputs:
            % D: DirectSolver
            % glob: GlobalOutside
            % patches: Patches
            % interfaces: Interfaces
            %
            % Outputs:
            % If change of variable,
            % u0: (mU+m_w+m_l)-by-1 doubles containing initial multiscale solution u0=(U0,w0,lambda0)
            % v0: (mU+m_w+m_l)-by-1 doubles containing initial multiscale velocity v0=(vU0,vw0,vlambda0)
            % else,
            % u0: (mU+m_z)-by-1 doubles containing initial multiscale solution u0=(U0,z0)
            % v0: (mU+m_z)-by-1 doubles containing initial multiscale velocity v0=(vU0,vz0)
            % where
            % n is the number of patches
            % m_U is the dimension of the spatial approximation space of global solution U
            % m_w is the dimension of the spatial approximation space of local solution w
            % m_l is the dimension of the spatial approximation space of Lagrange multiplier lambda
            % m_z is the dimension of the spatial approximation space of local solution z
            
            n = numel(patches);
            
            patch = patches.patches;
            interface = interfaces.interfaces;
            
            sz_U = getnbddlfree(glob.S);
            if ~D.changeOfVariable
                % Without change of variable
                sz_w = cellfun(@(patch) getnbddl(patch.S),patch,'UniformOutput',false);
                sz_lambda = cellfun(@(interface) getnbddl(interface.S),interface,'UniformOutput',false);
                % Multiscale initial solution vector u0=(U0,w0{k},lambda0{k})
                u0 = cell(1+2*n,1);
                if isempty(glob.u0)
                    u0{1} = zeros(sz_U,1);
                else
                    u0{1} = glob.u0;
                end
                for k=1:n
                    if isempty(patch{k}.u0)
                        u0{2*k} = zeros(sz_w{k},1);
                    else
                        u0{2*k} = patch{k}.u0;
                    end
                    if isempty(interface{k}.u0)
                        u0{2*k+1} = zeros(sz_lambda{k},1);
                    else
                        u0{2*k+1} = interface{k}.u0;
                    end
                end
                
                if nargout>1
                    % Multiscale initial velocity vector v0=(vU0,vw0{k},vlambda0{k})
                    v0 = cell(1+2*n,1);
                    if isempty(glob.v0)
                        v0{1} = zeros(sz_U,1);
                    else
                        v0{1} = glob.v0;
                    end
                    for k=1:n
                        if isempty(patch{k}.v0)
                            v0{2*k} = zeros(sz_w{k},1);
                        else
                            v0{2*k} = patch{k}.v0;
                        end
                        if isempty(interface{k}.v0)
                            v0{2*k+1} = zeros(sz_lambda{k},1);
                        else
                            v0{2*k+1} = interface{k}.v0;
                        end
                    end
                end
            else
                % With change of variable
                sz_z = cellfun(@(patch) getnbddlfree(patch.S),patch,'UniformOutput',false);
                % Multiscale initial solution vector u0=(U0,z0{k})
                u0 = cell(1+n,1);
                if isempty(glob.u0)
                    u0{1} = zeros(sz_U,1);
                else
                    u0{1} = glob.u0;
                end
                for k=1:n
                    if isempty(patch{k}.u0)
                        u0{k+1} = zeros(sz_z{k},1);
                    else
                        u0{k+1} = freevector(patch{k}.S,patch{k}.u0);
                    end
                end
                
                if nargout>1
                    % Multiscale initial solution vector v0=(vU0,vz0{k})
                    v0 = cell(1+n,1);
                    if isempty(glob.v0)
                        v0{1} = zeros(sz_U,1);
                    else
                        v0{1} = glob.v0;
                    end
                    for k=1:n
                        if isempty(patch{k}.v0)
                            v0{k+1} = zeros(sz_z{k},1);
                        else
                            v0{k+1} = freevector(patch{k}.S,patch{k}.v0);
                        end
                    end
                end
            end
            u0 = cell2mat(u0);
            if nargout>1
                v0 = cell2mat(v0);
            end
        end
        
        function [u,U,w,lambda] = getSolution(D,glob,patches,interfaces,u)
            % function [u,U,w,lambda] = getSolution(D,glob,patches,interfaces,u)
            % function [u,U,w] = getSolution(D,glob,patches,interfaces,u)
            
            patch = patches.patches;
            interface = interfaces.interfaces;
            
            sz_U = getnbddlfree(glob.S);
            if ~D.changeOfVariable
                % Without change of variable
                sz_w = cellfun(@(patch) getnbddl(patch.S),patch,'UniformOutput',false);
                sz_lambda = cellfun(@(interface) getnbddl(interface.S),interface,'UniformOutput',false);
                rep = [sz_w{:};sz_lambda{:}];
                rep = [sz_U;rep(:)];
                % Multiscale solution u=(U,w{k},lambda{k})
                if isempty(D.timeSolver)
                    % Time-independent problem
                    u = mat2cell(u,rep);
                    U = u{1};
                    w = cellfun(@(patch) u{2*patch.number},patch,'UniformOutput',false);
                    lambda = cellfun(@(interface) u{2*interface.number+1},interface,'UniformOutput',false);
                else
                    % Time-dependent problem
                    T = gettimemodel(u);
                    u_val = getvalue(u);
                    u = mat2cell(u_val,rep);
                    U = TIMEMATRIX(u{1},T);
                    w = cellfun(@(patch) TIMEMATRIX(u{2*patch.number},T),patch,'UniformOutput',false);
                    lambda = cellfun(@(interface) TIMEMATRIX(u{2*interface.number+1},T),interface,'UniformOutput',false);
                end
            else
                % With change of variable
                sz_z = cellfun(@(patch) getnbddlfree(patch.S),patch,'UniformOutput',false);
                rep = [sz_z{:}];
                rep = [sz_U;rep(:)];
                % Multiscale solution u=(U,z{k})
                if isempty(D.timeSolver)
                    % Time-independent problem
                    u = mat2cell(u,rep);
                    U = u{1};
                    z = cellfun(@(patch) u{1+patch.number},patch,'UniformOutput',false);
                else
                    % Time-dependent problem
                    T = gettimemodel(u);
                    u_val = getvalue(u);
                    u = mat2cell(u_val,rep);
                    U = TIMEMATRIX(u{1},T);
                    z = cellfun(@(patch) TIMEMATRIX(u{1+patch.number},T),patch,'UniformOutput',false);
                end
                
                P = cellfun(@(interface) interface.P_patch'*interface.P_globOut,interface,'UniformOutput',false);
                % Local solution w{k}
                w = cellfun(@(patch,z,P) unfreevector(patch.S,z) + P*U,patch,z,P,'UniformOutput',false);
            end
        end
        
        function lambda = computeLagrangeMultiplier(D,patches,interfaces,w,vw,aw)
            % function lambda = computeLagrangeMultiplier(D,patches,interfaces,w)
            % function lambda = computeLagrangeMultiplier(D,patches,interfaces,w,vw)
            % function lambda = computeLagrangeMultiplier(D,patches,interfaces,w,vw,aw)
            
            patch = patches.patches;
            interface = interfaces.interfaces;
            
            if isempty(D.timeSolver)
                lambda = cellfun(@(patch,interface,w) computeLagrangeMultiplier(patch,interface,w),patch,interface,w,'UniformOutput',false);
            else
                if D.timeOrder==1
                    lambda = cellfun(@(patch,interface,w,vw) computeLagrangeMultiplier(patch,interface,w,vw),patch,interface,w,vw,'UniformOutput',false);
                elseif D.timeOrder==2
                    lambda = cellfun(@(patch,interface,w,vw,aw) computeLagrangeMultiplier(patch,interface,w,vw,aw),patch,interface,w,vw,aw,'UniformOutput',false);
                end
            end
        end
        
    end
    
end
