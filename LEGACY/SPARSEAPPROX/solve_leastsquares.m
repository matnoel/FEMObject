function varargout = solve_leastsquares(A,y,regul,param,cv,k,varargin)
% function [u,err] = solve_leastsquares(A,y,regul,param,cv,k,varargin)
% Solve (regularized) least-squares minimization problem
% and compute relative cross-validation error
% A : matrix containing the evaluations of basis functions
% y : response vector or matrix
% regul : type of regularization ('' or 'l0' or 'l1'), '' by default
% param : parameters for regularization, struct() by default
% cv : type of cross-validation procedure ('leaveout' or 'kfold'), 'leaveout' by default
% k : number of folds (only for k-fold cross-validation procedure), 10 by default
% 
% function [u,err] = solve_leastsquares(A,y)
% ordinary least-squares (OLS) minimization (without regularization)
% min_{ u } ||y-A*u||_2^2
% 
% function [u,err] = solve_leastsquares(A,y,'l0',param)
% least-squares minimization with l0-regularization (l0-pseudo-norm)
% min_{ u } ||y-A*u||_2^2 such that ||u||_0 <= L
% or
% min_{ u } ||u||_0 such that ||y-A*u||_2^2 <= eps
% or
% min_{ u } 0.5*||y-A*u||_2^2 + lambda*||u||_0
% parameters of the optimization procedure
% param.L : optional, maximum number of elements in each decomposition (not more than param.L non-zeros coefficients), min(N,P) by default
% param.eps : optional, threshold on the squared l2-norm of the residual, 0 by default
% param.lambda : optional, penalty parameter, 0 by default
% param.numThreads : optional, number of threads for exploiting  multi-core/multi-cpus, -1 by default, which automatically selects all the available cpus/cores
% Warning: All the columns of A should have unit norm !
% 
% function [u,err] = solve_leastsquares(A,y,'l1',param)
% least-squares minimization with l1-regularization (l1-norm)
% if param.mode = 0 : min_{ u } ||y-A*u||_2^2 such that ||u||_1 <= lambda
% if param.mode = 1 : min_{ u } ||u||_1 such that ||y-A*u||_2^2 <= lambda
% if param.mode = 2 : min_{ u } 0.5*||y-A*u||_2^2 + lambda*||u||_2 + 0.5*lambda2*||u||_2^2
% parameters of the optimization procedure
% param.mode : mode, 2 by default
% param.lambda : parameter selected by cross-validation technique, 0 by default
% param.lambda2 : optional parameter; for mode=0 and mode=1, it adds a ridge  on the Gram matrix
% param.numThreads : optional, number of threads for exploiting  multi-core/multi-cpus, -1 by default, which automatically selects all the available cpus/cores
% param.pos : optional, adds non-negativity constraints on the coefficients, false by default
% param.cholesky : optional, choose between Cholesky implementation or one based on the matrix inversion Lemma, false by default
% param.ols : optional, perform an orthogonal projection before returning the solution, false by default
% param.max_length_path : optional, maximum length of the path, 4*size(A,2) by default
%
% function [u,err] = solve_leastsquares(A,y,'l2',param)
% least-squares minimization with l2-regularization (l2-norm)
% min_{ u } 0.5*||y-A*u||_2^2 + lambda*||u||_2
% parameters of the optimization procedure
% param.lambda : regularization parameter
% param.numThreads : optional, number of threads for exploiting  multi-core/multi-cpus, -1 by default, which automatically selects all the available cpus/cores
% param.itermax : optional, maximum number of iterations, 100 by default
% param.tol : optional, tolerance for stopping criteration, which is a relative duality gap if it is available, or a relative change of parameters
%
% See also SPARSEAPPROX/solve_leastsquares_KL

if nargin<3 || isempty(regul)
    regul = '';
end
if nargin<4 || isempty(param)
    param = struct();
end
if nargin<5 || isempty(cv)
    cv = 'leaveout';
end
if strcmp(cv,'kfold') && (nargin<6 || isempty(k))
    k = 10;
end

N = size(y,1);
n = size(y,2);
P = size(A,2);

% Parameters of the optimization procedure
switch regul
    case 'l0'
        % param.L = min(N,P);
        % param.eps = 0.;
        if ~isfield(param,'lambda') || isempty(param.lambda)
            param.lambda = 0.;
        end
        % param.numThreads = -1;
    case 'l1'
        if ~isfield(param,'mode') || isempty(param.mode)
            param.mode = 2;
        end
        if ~isfield(param,'lambda') || isempty(param.lambda)
            param.lambda = 0.;
        end
        % param.lambda2 = 0.;
        % param.numThreads = -1;
        % param.pos = false;
        % param.cholesky = false;
        % param.ols = false;
        % param.max_length_path = 4*P;
    case 'l2'
        if ~isfield(param,'lambda') || isempty(param.lambda)
            param.lambda = 0.;
        end
        % param.numThreads = -1;
        % param.itermax = 100;
        % param.tol = 1e-4;
        
end

u = zeros(P,n);
err = zeros(1,n);
if isempty(regul)
    % if N<P
    %     error('The number of samples N = %d must be superior or equal to the number of unknown coefficients P = %d',N,P)
    % else
    % Resolution
    u = solve(A'*A,A'*y);
    % parfor i=1:n
    %     u(:,i) = regress(y(:,i),A);
    % end
    % Cross-validation error
    if nargout>1
        err = calc_fast_cv_error(y,A,u,cv,k,varargin{:});
%         fun = @(y,A) solve(A'*A,A'*y);
%         err = calc_cv_error(y,A,fun,cv,k,varargin{:});
    end
    % end
else
    if n==1
        switch regul
            case 'l0'
                % Warning: All the columns of A_OMP should have unit norm !
                D = sqrt(sum(A.^2));
                A_OMP = A*spdiags(1./D(:),0,length(D),length(D));
                % A_OMP = A./repmat(D,[N,1]);
                % Resolution using OMP algorithm
                [u(:,1),solpath] = mexOMP(full(y),full(A_OMP),param);
                u(:,1) = spdiags(D(:),0,length(D),length(D))*u(:,1);
                u(:,1) = full(u(:,1));
            case 'l1'
                % Resolution using LARS algorithm
                [u(:,1),solpath] = mexLasso(full(y),full(A),param);
            case 'l2'
                % Resolution using ridge regression solver with a conjugate gradient solver
                u(:,1) = mexRidgeRegression(full(y),full(A),zeros(size(A,2),1),param);
            otherwise
                error(['Regularization technique ' regul ' not implemented'])
        end
        if exist('solpath','var') && norm(solpath)~=0
            % Selection of optimal solution using cross-validation technique
            [u(:,1),err] = select_opt_path(y,A,solpath,cv,k,'parallel');
        end
    else
        parfor i=1:n % for each component
            switch regul
                case 'l0'
                    % Warning: All the columns of A_OMP should have unit norm !
                    D = sqrt(sum(A.^2));
                    A_OMP = A*spdiags(1./D(:),0,length(D),length(D));
                    % A_OMP = A./repmat(D,[N,1]);
                    % Resolution using OMP algorithm
                    [u(:,i),solpath] = mexOMP(full(y(:,i)),full(A_OMP),param);
                    u(:,i) = spdiags(D(:),0,length(D),length(D))*u(:,i);
                    u(:,i) = full(u(:,i));
                    if norm(solpath)~=0
                        % Selection of optimal solution using cross-validation technique
                        [u(:,i),err(i)] = select_opt_path(y(:,i),A,solpath,cv,k);
                    end
                case 'l1'
                    % Resolution using LARS algorithm
                    [u(:,i),solpath] = mexLasso(full(y(:,i)),full(A),param);
                    if norm(solpath)~=0
                        % Selection of optimal solution using cross-validation technique
                        [u(:,i),err(i)] = select_opt_path(y(:,i),A,solpath,cv,k);
                    end
                case 'l2'
                    % Resolution using ridge regression solver with a conjugate gradient solver
                    u(:,i) = mexRidgeRegression(full(y(:,i)),full(A),zeros(size(A,2),1),param);
                otherwise
                    error(['Regularization technique ' regul ' not implemented'])
            end
        end
    end
end

varargout{1} = u;
if nargout>1
    varargout{2} = err;
end

end

