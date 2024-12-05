function [u,err,x,PC_seq,err_seq,x_seq] = decompmatrix_leastsquares_KL(apc,N,fun,n,regul,param,cv,k,tolKL,cvKL,kKL,opts,varargin)
% function [u,err,x,PC_seq,err_seq,x_seq] = decompmatrix_leastsquares_KL(PC,N,fun,n,regul,param,cv,k,tolKL,cvKL,kKL,opts,varargin)
% Decomposition of a function of random variables on polynomial chaos basis
% by (regularized) least-squares minimization (or regression) using a priori empirical Karhunen-Loeve (KL) decomposition
% and cross-validation technique
% PC : (initial) POLYCHOAS or PCMODEL or FESTOMODEL or PCMATRIX or PCARRAY
% N : (initial) number of samples (regression points), length(PC) by default
% fun : function_handle of random variables RANDVARS(PC) associated to PC
% n : number of output arguments for fun
% regul : type of regularization ('' or 'l0' or 'l1'), '' by default
% param : parameters for regularization, struct() by default
% cv : type of cross-validation procedure ('leaveout' or 'kfold'), 'leaveout' by default
% k : number of folds (only for cv = 'kfold'), 10 by default
% tolKL : prescribed tolerance for truncated empirical KL decomposition
% cvKL : type of cross-validation procedure ('leaveout' or 'kfold'), 'leaveout' by default
% kKL : number of folds (only for k-fold cross-validation procedure), 10 by default
% 'ploterrorKL' in varargin : plot relative SVD truncation error w.r.t. SVD rank (true or false), false by default
% opts : struct of options for adaptive strategy
%        sampling : sampling strategy ('fixed' or 'adaptive'), 'fixed' by default
%        addsample : percentage of additional samples if 0 < addsample < 1 , 0.1 by default
%                    number of additional samples if addsample > 1
%        maxsample : maximal number of samples, Inf by default
%        basis : construction of PC basis ('fixed' or 'adaptive'), 'fixed' by default
%        algorithm : adaptive algorithm for the construction of a nested sequence of finite monotone/lower multi-index sets ('TP' or 'PD' or 'TD' or 'MS' or 'RMS'), 'RMS' by default
%                    'TP' or 'PD': isotropic Tensor Product (or Partial Degree) polynomial space
%                           multidimensional space of polynomials of partial degree less or equal to p (in each variable)
%                           update the partial degree by 1 in each dimension at each iteration
%                    'TD':  isotropic Total Degree polynomial space
%                           multidimensional space of polynomials of total degree less or equal to p
%                           update the total degree by 1 at each iteration
%                    'MS':  Margin Search strategy
%                           add a smallest monotone/lower subset S_n of the margin M_n of a given monotone/lower set A_n
%                           for which energy(S_n)>=bulkparam*energy(M_n), where bulkparam is a bulk parameter
%                    'RMS': Reduced Margin Search strategy
%                           add a smallest monotone/lower subset S_n of the reduced margin M_n of a given monotone/lower set A_n
%                           for which energy(S_n)>=bulkparam*energy(M_n), where bulkparam is a bulk parameter
%        bulkparam : bulk parameter in (0,1) such that energy(S_n)>=bulkparam*energy(M_n), 0.5 by default
%                    bulkparam = 1 selects all multi-indices in the (reduced) margin M_n
%                    bulkparam = 0 selects the multi-index in the (reduced) margin M_n
%                                  corresponding to the maximal norm of the expansion coefficients
%        tol : prescribed tolerance for cross-validation error, 1e-4 by default
%        tolstagn : prescribed stagnation tolerance for cross-validation error, 5e-2 by default
%        toloverfit : prescribed tolerance to detect overfitting for cross-validation error such that err>=toloverfit*err_old, 1.1 by default
%        correction : correction for cross-validation error, true by default
% varargin : arguments of function fun
% call y = fun(x,varargin{:}) for matrix sample x
% x(:,i) is the sample vector containing all the evaluations of random variable i
% x(k,:) is the sample vector containing the k-th evaluation of all random variables
% y is the response matrix contaning the values returned by function fun applied to sample matrix x
% if PC is a POLYCHAOS, x is the sample matrix containing the evaluations of random variables associated to POLYCHAOS
%          a PCMODEL or a FESTOMODEL, x is the sample matrix containing the evaluations of random variables associated to PCMODEL
%          a PCMATRIX or a PCARRAY, x is the sample matrix containing the evaluations of random variables associated to PCMATRIX
% 
% function [u,err,x,PC_seq,err_seq,x_seq] = decompmatrix_leastsquares_KL(PC,N,fun)
% ordinary least-squares (OLS) minimization (without regularization)
% min_{ u } ||y-A*u||_2^2
% 
% function [u,err,x,PC_seq,err_seq,x_seq] = decompmatrix_leastsquares_KL(PC,N,fun,'l0',param)
% least-squares minimization with l0-regularization (l0-pseudo-norm)
% min_{ u } ||y-A*u||_2^2 such that ||u||_0 <= L
% or
% min_{ u } ||u||_0 such that ||y-A*u||_2^2 <= eps
% or
% min_{ u } 0.5*||y-A*u||_2^2 + lambda*||u||_0
% param : struct of parameters of the optimization procedure
%         L : optional, maximum number of elements in each decomposition (not more than param.L non-zeros coefficients), min(N,P) by default
%         eps : optional, threshold on the squared l2-norm of the residual, 0 by default
%         lambda : optional, penalty parameter, 0 by default
%         numThreads : optional, number of threads for exploiting  multi-core/multi-cpus, -1 by default, which automatically selects all the available cpus/cores of the machine
% Warning: All the columns of A should have unit norm !
% 
% function [u,err,x,PC_seq,err_seq,x_seq] = decompmatrix_leastsquares_KL(PC,N,fun,'l1',param)
% least-squares minimization with l1-regularization (l1-norm)
% if param.mode = 0 : min_{ u } ||y-A*u||_2^2 such that ||u||_1 <= lambda
% if param.mode = 1 : min_{ u } ||u||_1 such that ||y-A*u||_2^2 <= lambda
% if param.mode = 2 : min_{ u } 0.5*||y-A*u||_2^2 + lambda*||u||_1 + 0.5*lambda2*||u||_2^2
% param : struct of parameters of the optimization procedure
%         mode : mode, 2 by default
%         lambda : parameter selected by cross-validation technique, 0 by default
%         lambda2 : optional parameter; for mode=0 and mode=1, it adds a ridge  on the Gram matrix
%         numThreads : optional, number of threads for exploiting  multi-core/multi-cpus, -1 by default, which automatically selects all the available cpus/cores of the machine
%         pos : optional, adds non-negativity constraints on the coefficients, false by default
%         cholesky : optional, choose between Cholesky implementation or one based on the matrix inversion Lemma, false by default
%         ols : optional, perform an orthogonal projection before returning the solution, false by default
%         max_length_path : optional, maximum length of the path, 4*size(A,2) by default
%
% function [u,err,x,PC_seq,err_seq,x_seq] = decompmatrix_leastsquares(PC,N,fun,'l2',param)
% least-squares minimization with l2-regularization (l2-norm)
% min_{ u } 0.5*||y-A*u||_2^2 + lambda*||u||_2
% parameters of the optimization procedure
% param.lambda : regularization parameter
% param.numThreads : optional, number of threads for exploiting  multi-core/multi-cpus, -1 by default, which automatically selects all the available cpus/cores
% param.itermax : optional, maximum number of iterations, 100 by default
% param.tol : optional, tolerance for stopping criteration, which is a relative duality gap if it is available, or a relative change of parameters

if nargin<4 || isempty(n)
    n = 1;
end
if nargin<5 || isempty(regul)
    regul = '';
end
if nargin<6 || isempty(param)
    param = struct();
end
if nargin<7 || isempty(cv)
    cv = 'leaveout';
end
if nargin<8 || isempty(k)
    k = 10;
end
if nargin<9 || isempty(tolKL)
    tolKL = [];
end
if nargin<10 || isempty(cvKL)
    cvKL = 'leaveout';
end
if nargin<11 || isempty(kKL)
    kKL = 10;
end
if nargin<12 || ~isfield(opts,'sampling')
    opts.sampling = 'fixed';
end
if nargin<12 || ~isfield(opts,'addsample')
    opts.addsample = 0.1;
end
if nargin<12 || ~isfield(opts,'maxsample')
    opts.maxsample = Inf;
end
if nargin<12 || ~isfield(opts,'basis')
    opts.basis = 'fixed';
end
if nargin<12 || ~isfield(opts,'maxcoeff')
    opts.maxcoeff = Inf;
end
if nargin<12 || ~isfield(opts,'algorithm')
    opts.algorithm = 'RMS';
end
if nargin<12 || ~isfield(opts,'bulkparam')
    opts.bulkparam = 0.5;
end
if nargin<12 || ~isfield(opts,'tol')
    opts.tol = 1e-4;
end
if nargin<12 || ~isfield(opts,'tolstagn')
    opts.tolstagn = 5e-2;
end
if nargin<12 || ~isfield(opts,'toloverfit')
    opts.toloverfit = 1.1;
end
if nargin<12 || ~isfield(opts,'correction')
    opts.correction = true;
end

if ischarin('display',varargin)
    display_ = 1;
    varargin = delonlycharin('display',varargin);
else
    display_ = 0;
end
if ischarin('displayiter',varargin)
    display_iter = 1;
    varargin = delonlycharin('displayiter',varargin);
else
    display_iter = 0;
end

fun = fcnchk(fun);
PC = getPC(apc);

if isempty(N)
    N = length(PC);
elseif isa(N,'function_handle')
    N = N(PC);
end
P = length(PC);
if isempty(regul) && N<P
    warning('The number of samples N = %d must be superior or equal to the number of unknown coefficients P = %d. N will be set to the value of P.',N,P);
    N = P;
end

fP = '%5d';
fN = '%5d';
fr = '%8d';
ferr = '%14.4e';
if display_ || display_iter
    fprintf('\n+-------+-------+----------+----------------+\n');
    fprintf('|   P   |   N   | SVD rank | Relative error |\n');
    fprintf('+-------+-------+----------+----------------+\n');
end

[A,xpc] = random(PC,N);
x = transfer(apc,xpc);
[y,s] = calc_response(fun,x,n,varargin{:});

u = cell(1,n);
err = cell(1,n);
r = cell(1,n);
parfor i=1:n
    if isa(y{i},'cell')
        for l=1:length(y{i})
            if opts.correction
                [u{i}{l},err{i}{l},r{i}{l}] = solve_leastsquares_KL(A,y{i}{l}',regul,param,cv,k,tolKL,cvKL,kKL,'correction');
            else
                [u{i}{l},err{i}{l},r{i}{l}] = solve_leastsquares_KL(A,y{i}{l}',regul,param,cv,k,tolKL,cvKL,kKL);
            end
            u{i}{l} = u{i}{l}';
            if display_iter
                fprintf(['| ' fP ' | ' fN ' | ' fr ' | ' ferr ' |\n'],P,N,r{i}{l},norm(err{i}{l}));
            end
        end
    else
        if opts.correction
            [u{i},err{i},r{i}] = solve_leastsquares_KL(A,y{i}',regul,param,cv,k,tolKL,cvKL,kKL,'correction');
        else
            [u{i},err{i},r{i}] = solve_leastsquares_KL(A,y{i}',regul,param,cv,k,tolKL,cvKL,kKL);
        end
        u{i} = u{i}';
        if display_iter
            fprintf(['| ' fP ' | ' fN ' | ' fr ' | ' ferr ' |\n'],P,N,r{i},norm(err{i}));
        end
    end
end
if display_iter && ~(n==1 && ~isa(y{1},'cell'))
    fprintf( '|       |       |          |                |\n');
end

PC_seq = initialize_seq(PC,y);
err_seq = initialize_seq(err,y);
x_seq{1} = x;

sampling = opts.sampling;
if strcmp(sampling,'adaptive')
    err_old = initialize(Inf,y);
    while ~convergence(err,opts.tol) && ~stagnation(err,err_old,opts.tolstagn)
        if N>=opts.maxsample
            sampling = 'fixed';
            break
        end
        if opts.addsample<1
            N_add = ceil(opts.addsample*N);
        else
            N_add = opts.addsample;
        end
        N = N + N_add;
        [A_add,xpc_add] = random(PC,N_add);
        xpc = [xpc;xpc_add];
        A = [A;A_add];
        x_add = transfer(apc,xpc_add);
        x = [x;x_add];
        x_seq{end+1} = x;
        y_add = calc_response(fun,x_add,n,varargin{:});
        err_old = err;
        for i=1:n
            if isa(y{i},'cell')
                for l=1:length(y{i})
                    y{i}{l} = [y{i}{l} y_add{i}{l}];
                    if opts.correction
                        [u{i}{l},err{i}{l},r{i}{l}] = solve_leastsquares_KL(A,y{i}{l}',regul,param,cv,k,tolKL,cvKL,kKL,'correction');
                    else
                        [u{i}{l},err{i}{l},r{i}{l}] = solve_leastsquares_KL(A,y{i}{l}',regul,param,cv,k,tolKL,cvKL,kKL);
                    end
                    u{i}{l} = u{i}{l}';
                    PC_seq{i}{l}{end+1} = PC;
                    err_seq{i}{l}{end+1} = norm(err{i}{l});
                    if display_iter
                        fprintf(['|       | ' fN ' | ' fr ' | ' ferr ' |\n'],N,r{i}{l},norm(err{i}{l}));
                    end
                end
            else
                y{i} = [y{i} y_add{i}];
                if opts.correction
                    [u{i},err{i},r{i}] = solve_leastsquares_KL(A,y{i}',regul,param,cv,k,tolKL,cvKL,kKL,'correction');
                else
                    [u{i},err{i},r{i}] = solve_leastsquares_KL(A,y{i}',regul,param,cv,k,tolKL,cvKL,kKL);
                end
                u{i} = u{i}';
                PC_seq{i}{end+1} = PC;
                err_seq{i}{end+1} = norm(err{i});
                if display_iter
                    fprintf(['|       | ' fN ' | ' fr ' | ' ferr ' |\n'],N,r{i},norm(err{i}));
                end
            end
        end
        if display_iter && ~(n==1 && ~isa(y{1},'cell'))
            fprintf( '|       |       |          |                |\n');
        end
    end
end

basis = opts.basis;
if strcmp(basis,'adaptive')
    if ~islower(PC)
        error('The initial multi-index set is not a lower/downward closed/monotone set.')
    end
    rv = RANDVARS(PC);
    PC = initialize(PC,y);
    A = initialize(A,y);
    err = initialize(err);
    while ~convergence(err,opts.tol)
        if strcmp(basis,'fixed') && strcmp(sampling,'fixed')
            break
        end
        if strcmp(basis,'adaptive')
            err_old = initialize(Inf,y);
            while ~convergence(err,opts.tol)
                if maxcoeff(PC,opts.maxcoeff) || (strcmp(sampling,'fixed') && ~update_basis(PC,err,err_old,regul,N,opts))
                    basis = 'fixed';
                    break
                end
                if ~update_basis(PC,err,err_old,regul,N,opts)
                    break
                end
                x_seq{end+1} = x;
                for i=1:n
                    if isa(y{i},'cell')
                        for l=1:length(y{i})
                            if update_basis(PC{i}{l},err{i}{l},err_old{i}{l},regul,N,opts)
                                if opts.correction
                                    funsolve = @(A) solve_leastsquares_KL(A,y{i}{l}',regul,param,cv,k,tolKL,cvKL,kKL,'correction');
                                else
                                    funsolve = @(A) solve_leastsquares_KL(A,y{i}{l}',regul,param,cv,k,tolKL,cvKL,kKL);
                                end
                                PC{i}{l} = adapt(PC{i}{l},funsolve,xpc,opts.algorithm,opts.bulkparam);
                                A{i}{l} = polyval(PC{i}{l},xpc);
                                err_old{i}{l} = err{i}{l};
                                if opts.correction
                                    [u{i}{l},err{i}{l},r{i}{l}] = solve_leastsquares_KL(A{i}{l},y{i}{l}',regul,param,cv,k,tolKL,cvKL,kKL,'correction');
                                else
                                    [u{i}{l},err{i}{l},r{i}{l}] = solve_leastsquares_KL(A{i}{l},y{i}{l}',regul,param,cv,k,tolKL,cvKL,kKL);
                                end
                                u{i}{l} = u{i}{l}';
                                PC_seq{i}{l}{end+1} = PC{i}{l};
                                err_seq{i}{l}{end+1} = norm(err{i}{l});
                            end
                            if display_iter
                                P = length(PC{i}{l});
                                fprintf(['| ' fP ' |       | ' fr ' | ' ferr ' |\n'],P,r{i}{l},norm(err{i}{l}));
                            end
                        end
                    else
                        if update_basis(PC{i},err{i},err_old{i},regul,N,opts)
                            if opts.correction
                                funsolve = @(A) solve_leastsquares_KL(A,y{i}',regul,param,cv,k,tolKL,cvKL,kKL,'correction');
                            else
                                funsolve = @(A) solve_leastsquares_KL(A,y{i}',regul,param,cv,k,tolKL,cvKL,kKL);
                            end
                            PC{i} = adapt(PC{i},funsolve,xpc,opts.algorithm,opts.bulkparam);
                            A{i} = polyval(PC{i},xpc);
                            err_old{i} = err{i};
                            if opts.correction
                                [u{i},err{i},r{i}] = solve_leastsquares_KL(A{i},y{i}',regul,param,cv,k,tolKL,cvKL,kKL,'correction');
                            else
                                [u{i},err{i},r{i}] = solve_leastsquares_KL(A{i},y{i}',regul,param,cv,k,tolKL,cvKL,kKL);
                            end
                            u{i} = u{i}';
                            PC_seq{i}{end+1} = PC{i};
                            err_seq{i}{end+1} = norm(err{i});
                        end
                        if display_iter
                            P = length(PC{i});
                            fprintf(['| ' fP ' |       | ' fr ' | ' ferr ' |\n'],P,r{i},norm(err{i}));
                        end
                    end
                end
                if display_iter && ~(n==1 && ~isa(y{1},'cell'))
                    fprintf( '|       |       |          |                |\n');
                end
            end
        end
        if strcmp(sampling,'adaptive')
            err_old = initialize(Inf,y);
            while ~convergence(err,opts.tol)
                if N>=opts.maxsample || (strcmp(basis,'fixed') && stagnation(err,err_old,opts.tolstagn))
                    sampling = 'fixed';
                    break
                end
                if stagnation(err,err_old,opts.tolstagn)
                    break
                end
                if opts.addsample<1
                    N_add = ceil(opts.addsample*N);
                else
                    N_add = opts.addsample;
                end
                N = N + N_add;
                xpc_add = random(rv,N_add,1);
                xpc_add = [xpc_add{:}];
                xpc = [xpc;xpc_add];
                x_add = transfer(apc,xpc_add);
                x = [x;x_add];
                x_seq{end+1} = x;
                y_add = calc_response(fun,x_add,n,varargin{:});
                err_old = err;
                for i=1:n
                    if isa(y{i},'cell')
                        for l=1:length(y{i})
                            A_add = polyval(PC{i}{l},xpc_add);
                            A{i}{l} = [A{i}{l};A_add];
                            y{i}{l} = [y{i}{l} y_add{i}{l}];
                            if opts.correction
                                [u{i}{l},err{i}{l},r{i}{l}] = solve_leastsquares_KL(A{i}{l},y{i}{l}',regul,param,cv,k,tolKL,cvKL,kKL,'correction');
                            else
                                [u{i}{l},err{i}{l},r{i}{l}] = solve_leastsquares_KL(A{i}{l},y{i}{l}',regul,param,cv,k,tolKL,cvKL,kKL);
                            end
                            u{i}{l} = u{i}{l}';
                            PC_seq{i}{l}{end+1} = PC{i}{l};
                            err_seq{i}{l}{end+1} = norm(err{i}{l});
                            if display_iter
                                fprintf(['|       | ' fN ' | ' fr ' | ' ferr ' |\n'],N,r{i}{l},norm(err{i}{l}));
                            end
                        end
                    else
                        A_add = polyval(PC{i},xpc_add);
                        A{i} = [A{i};A_add];
                        y{i} = [y{i} y_add{i}];
                        if opts.correction
                            [u{i},err{i},r{i}] = solve_leastsquares_KL(A{i},y{i}',regul,param,cv,k,tolKL,cvKL,kKL,'correction');
                        else
                            [u{i},err{i},r{i}] = solve_leastsquares_KL(A{i},y{i}',regul,param,cv,k,tolKL,cvKL,kKL);
                        end
                        u{i} = u{i}';
                        PC_seq{i}{end+1} = PC{i};
                        err_seq{i}{end+1} = norm(err{i});
                        if display_iter
                            fprintf(['|       | ' fN ' | ' fr ' | ' ferr ' |\n'],N,r{i},norm(err{i}));
                        end
                    end
                end
                if display_iter && ~(n==1 && ~isa(y{1},'cell'))
                    fprintf( '|       |       |          |                |\n');
                end
            end
        end
    end
end

if display_ && ~display_iter
    for i=1:n
        if isa(err{i},'cell')
            for l=1:length(err{i})
                if strcmp(opts.basis,'adaptive')
                    P = length(PC{i}{l});
                else
                    P = length(PC);
                end
                fprintf(['| ' fP ' | ' fN ' | ' ferr ' |\n'],P,N,err{i}{l});
            end
        else
            if strcmp(opts.basis,'adaptive')
                P = length(PC{i});
            else
                P = length(PC);
            end
            fprintf(['| ' fP ' | ' fN ' | ' ferr ' |\n'],P,N,err{i});
        end
    end
end

if display_ || display_iter
    fprintf( '+-------+-------+----------+----------------+\n');
end

for i=1:n
    if isa(u{i},'cell')
        for l=1:length(u{i})
            if strcmp(opts.basis,'adaptive')
                u{i}{l} = PCMATRIX(u{i}{l},s{i}{l},PC{i}{l});
            else
                u{i}{l} = PCMATRIX(u{i}{l},s{i}{l},PC);
            end
        end
    else
        if strcmp(opts.basis,'adaptive')
            u{i} = PCMATRIX(u{i},s{i},PC{i});
        else
            u{i} = PCMATRIX(u{i},s{i},PC);
        end
    end
end

if n==1 && ~isa(u{1},'cell')
    u = [u{:}];
    err = [err{:}];
    PC_seq = [PC_seq{:}];
    err_seq = [err_seq{:}];
end

end

function a = initialize_seq(a_init,y)

if nargin==1
    n = length(a_init);
    a = cell(1,n);
    for i=1:n
        if isa(a_init{i},'cell')
            for l=1:length(a_init{i})
                a{i}{l}{1} = a_init{i}{l};
            end
        else
            a{i}{1} = a_init{i};
        end
    end
else
    n = length(y);
    a = cell(1,n);
    for i=1:n
        if isa(y{i},'cell')
            for l=1:length(y{i})
                a{i}{l}{1} = a_init;
            end
        else
            a{i}{1} = a_init;
        end
    end
end

end

function a = initialize(a_init,y)

if nargin==1
    a = a_init;
else
    n = length(y);
    a = cell(1,n);
    for i=1:n
        if isa(y{i},'cell')
            for l=1:length(y{i})
                a{i}{l} = a_init;
            end
        else
            a{i} = a_init;
        end
    end
end

end

function bool = finite(err)

bool = true;
if ~isa(err,'cell')
    if ~isfinite(norm(err))
        bool = false;
    end
else
    n = length(err);
    for i=1:n
        if isa(err{i},'cell')
            for l=1:length(err{i})
                if ~isfinite(norm(err{i}{l}))
                    bool = false;
                    return
                end
            end
        else
            if ~isfinite(norm(err{i}))
                bool = false;
                return
            end
        end
    end
end

end

function bool = convergence(err,tol)

bool = true;
if ~isa(err,'cell')
    if norm(err)>=tol
        bool = false;
    end
else
    n = length(err);
    for i=1:n
        if isa(err{i},'cell')
            for l=1:length(err{i})
                if norm(err{i}{l})>=tol
                    bool = false;
                    return
                end
            end
        else
            if norm(err{i})>=tol
                bool = false;
                return
            end
        end
    end
end

end

function bool = stagnation(err,err_old,tolstagn)

bool = true;
if ~isa(err,'cell')
    if norm(err-err_old)/norm(err)>=tolstagn || ~isfinite(norm(err))
        bool = false;
    end
else
    n = length(err);
    for i=1:n
        if isa(err{i},'cell')
            for l=1:length(err{i})
                if norm(err{i}{l}-err_old{i}{l})/norm(err{i}{l})>=tolstagn || ~isfinite(norm(err{i}{l}))
                    bool = false;
                    return
                end
            end
        else
            if norm(err{i}-err_old{i})/norm(err{i})>=tolstagn || ~isfinite(norm(err{i}))
                bool = false;
                return
            end
        end
    end
end

end

function bool = overfitting(err,err_old,toloverfit)

bool = false;
if ~isa(err,'cell')
    if norm(err)>=toloverfit*norm(err_old)
        bool = true;
    end
else
    n = length(err);
    for i=1:n
        if isa(err{i},'cell')
            for l=1:length(err{i})
                if norm(err{i}{l})>=toloverfit*norm(err_old{i}{l})
                    bool = true;
                    return
                end
            end
        else
            if norm(err{i})>=toloverfit*norm(err_old{i})
                bool = true;
                return
            end
        end
    end
end

end

function bool = maxcoeff(PC,P_max)

bool = true;
if ~isa(PC,'cell')
    P = length(PC);
    if P<P_max
        bool = false;
    end
else
    n = length(PC);
    for i=1:n
        if isa(PC{i},'cell')
            for l=1:length(PC{i})
                P = length(PC{i}{l});
                if P<P_max
                    bool = false;
                    return
                end
            end
        else
            P = length(PC{i});
            if P<P_max
                bool = false;
                return
            end
        end
    end
end

end

function bool = adapt_basis(PC,regul,N,algorithm)

if isempty(regul)
    bool = false;
    if ~isa(PC,'cell')
        switch algorithm
            case {'TP','PD','TD'}
                h = RANDPOLYS(PC);
                p = getorder(PC);
                if strcmp(algorithm,'TP') || strcmp(algorithm,'PD')
                    PC_new = POLYCHAOS(h,p+1,'typebase',2);
                elseif strcmp(algorithm,'TD')
                    PC_new = POLYCHAOS(h,p+1,'typebase',1);
                end
            case {'MS','RMS'}
                if strcmp(algorithm,'MS')
                    ind_marg = getindices(PC,'margin');
                elseif strcmp(algorithm,'RMS')
                    ind_marg = getindices(PC,'reducedmargin');
                end
                PC_new = addindices(PC,ind_marg,'update');
        end
        P_new = length(PC_new);
        if P_new<=N
            bool = true;
        end
    else
        n = length(PC);
        for i=1:n
            if isa(PC{i},'cell')
                for l=1:length(PC{i})
                    bool = adapt_basis(PC{i}{l},regul,N,algorithm);
                    if bool
                        return
                    end
                end
            else
                bool = adapt_basis(PC{i},regul,N,algorithm);
                if bool
                    return
                end
            end
        end
    end
else
    bool = true;
end

end

function bool = update_basis(PC,err,err_old,regul,N,opts)

bool = false;
if ~isa(PC,'cell')
    if ~maxcoeff(PC,opts.maxcoeff) && ~convergence(err,opts.tol) && ~overfitting(err,err_old,opts.toloverfit) && adapt_basis(PC,regul,N,opts.algorithm)% && ~stagnation(err,err_old,opts.tolstagn)
        bool = true;
    end
else
    n = length(PC);
    for i=1:n
        if isa(PC{i},'cell')
            for l=1:length(PC{i})
                bool = update_basis(PC{i}{l},err{i}{l},err_old{i}{l},regul,N,opts);
                if bool
                    return
                end
            end
        else
            bool = update_basis(PC{i},err{i},err_old{i},regul,N,opts);
            if bool
                return
            end
        end
    end
end

end
