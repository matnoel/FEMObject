function [sol,modelerror]=sparse_least_square_solution(param,fs,rs,fsModel,rsModel,RV)
% function [sol,modelerror]=sparse_least_square_solution(param,fs,rs,fsModel,rsModel,RV)

p = param.order;
if (strcmp(param.basistype,'total'))
    PC = POLYCHAOS(RV,p,'pcg','typebase',1); % PC basis
elseif (strcmp(param.basistype,'partial'))
    PC = POLYCHAOS(RV,p,'pcg','typebase',2);
end
%rs_pc = transfer(RV,RANDVARS(PC),rs); % mapping of each random variable RV (defined over [-pi,pi])
                                      % into RANDVARS(PC) (defined over [-1,1] for Legendre polynomials)
                                      % for each random sample
A = polyval(PC,rs); % evaluation of each polynomial of PC basis on random samples (N-by-P matrix)

% Type of regularization
% 'l0' : regularization with respect to l0-norm
% 'l1' : regularization with respect to l1-norm
% '' : no regularization
regul = 'l1';

N = size(D,1);
P = size(D,2);
if strcmp(regul,'l0')
    % Least-squares minimization problem with l0-regularization
    % min_{ x } ||y-A*x||_2^2 s.t. ||x||_0 <= L
    % or
    % min_{ x } ||x||_0 s.t. ||y-A*x||_2^2 <= eps
    % or
    % min_{ x } 0.5*||y-A*x||_2^2 + lambda*||x||_0
    
    % Parameters of the optimization procedure
    % param.L = min(N,P); % maximum number of elements in each decomposition (not more than param.L non-zeros coefficients);
    % the default choice is min(N,P)
    % param.eps = 0.; % threshold on the squared l2-norm of the residual
    % the default choice is 0
    param.lambda = 0.; % penalty parameter;
    % the default choice is 0
    % param.numThreads = -1; % number of threads for exploiting  multi-core/multi-cpus;
    % the default choice is -1, which automatically selects all the available cpus/cores of the machine
    
    % Resolution using OMP algorithm
    A_OMP = A./repmat(sqrt(sum(D.^2)), [N,1]); % Warning: All the columns of A_OMP should have unit norm !
    [x,path] = mexLasso(fs,full(A_OMP),param);
    % Cross-validation technique
    [path_opt,epsil] = select_opt_path(full(A),fs,path);
    % Solution
    sol = PCMATRIX(path_opt,[1,1],PC); % P-by-1 vector
    
elseif strcmp(regul,'l1')
    % Least-squares minimization problem with l1-regularization
    % when param.mode = 0 : min_{ x } ||y-A*x||_2^2 s.t. ||x||_1 <= lambda
    % when param.mode = 1 : min_{ x } ||x||_1 s.t. ||y-A*x||_2^2 <= lambda
    % when param.mode = 2 : min_{ x } 0.5*||y-A*x||_2^2 + lambda*||x||_1 + 0.5*lambda2*||x||_2^2
    
    % Parameters of the optimization procedure
    param.mode = 2; % the default choice is 2
    param.lambda = 0.; % parameter
    % param.lambda2 = 0.; % optional parameter; for mode=0 and mode=1, it adds a ridge  on the Gram matrix
    % param.numThreads = -1; % number of threads for exploiting  multi-core/multi-cpus;
    % the default choice is -1, which automatically selects all the available cpus/cores of the machine
    
    % Resolution using LARS algorithm
    [path_opt,path] = mexLasso(fs,full(A),param);
    % Cross-validation technique
    [path_opt,epsil] = select_opt_path(full(A),fs,path);
    % Solution
    sol = PCMATRIX(path_opt,[1,1],PC);
    
else
    % Least-squares minimization problem without regularization
    % min_{ x } ||y-A*x||_2^2
    
    % Resolution
    x = solve(A'*A,A'*fs);
    % Solution
    sol = PCMATRIX(x,[1,1],PC);
end
error = norm(randomeval(sol,rsModel)'-fsModel)/norm(fsModel);
modelerror = error;
