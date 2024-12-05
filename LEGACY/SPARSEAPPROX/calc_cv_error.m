function [err,delta] = calc_cv_error(y,A,fun,cv,k,varargin)
% function [err,delta] = calc_cv_error(y,A,fun,cv,k)
% Compute relative cross-validation error for coefficients vector/matrix u
% = fun(y,A) and residuals delta
% A : matrix containing the evaluations of basis functions
% y : response vector or matrix
% fun : function handle
% cv : type of cross-validation procedure ('leaveout' or 'kfold'), 'leaveout' by default
% k : number of folds (only for k-fold cross-validation procedure), min(10,N) by default where N is the number of samples
% 
% function [err,delta] = calc_cv_error(y,A,fun,'leaveout')
% Compute relative leave-one-out cross-validation error for coefficients
% vector/matrix u = fun(y,A) and residuals delta
% 
% function [err,delta] = calc_cv_error(y,A,fun,'kfold',k)
% Compute relative k-fold cross-validation error for coefficients
% vector/matrix u = fun(y,A) and residuals delta
% 
% function [err,delta] = calc_cv_error(y,A,fun,'leaveout',[],'correction')
% Compute corrected relative leave-one-out cross-validation error for
% coefficients vector/matrix u = fun(y,A) and residuals delta
% 
% function [err,delta] = calc_cv_error(y,A,fun,'kfold',k,'correction')
% Compute corrected relative k-fold cross-validation error for coefficients
% vector/matrix u = fun(y,A) and residuals delta
% 
% See also SPARSEAPPROX/calc_fast_cv_error

N = size(y,1);
n = size(y,2);
P = size(A,2);

if nargin<4 || isempty(cv)
    cv = 'leaveout';
end
if strcmp(cv,'kfold') && (nargin<5 || isempty(k))
    k = min(10,N);
end

% initialization
err = zeros(1,n);
delta = zeros(size(y));

% create a random partition of nearly equal size for leave-one-out or k-fold cross-validation on N observations
switch cv
    case 'leaveout'
        cvp = cvpartition(N,'leaveout');
    case 'kfold'
        cvp = cvpartition(N,'kfold',min(k,N));
    otherwise
        error(['Cross-validation ' cv ' not implemented'])
end

% check whether there are enough samples when performing ordinary least-squares minimization
f = functions(fun);
regul = f.workspace{2}.regul;
if isempty(regul) && any(cvp.TrainSize<P)
    % warning('Not egough samples for performing OLS on the training set')
    err(:) = Inf;
    delta(:) = Inf;
    return
end

% compute absolute cross-validation error (cross-validation error estimate of mean-squared error)
% also called mean predicted residual sum of squares (PRESS) or empirical mean-squared predicted residual
errors = cell(cvp.NumTestSets,1);
deltas = cell(cvp.NumTestSets,1);
parfor i=1:cvp.NumTestSets % for each fold i
    trIdx = cvp.training(i);
    teIdx = cvp.test(i);
    Atrain = A(trIdx,:);
    Atest = A(teIdx,:);
    ytrain = y(trIdx,:);
    ytest = y(teIdx,:);
    % compute predicted solution for training set i
    utrain = fun(ytrain,Atrain);
    % compute predicted residual over test set i
    deltas{i} = ytest-Atest*utrain;
    % compute absolute cross-validation error for fold i
    errors{i} = sum(deltas{i}.^2,1)/cvp.TestSize(i);
end
for i=1:cvp.NumTestSets
    teIdx = cvp.test(i);
    delta(teIdx,:) = deltas{i};
end
% average over k = cvp.NumTestSets folds
err = sum(cell2mat(errors),1)/cvp.NumTestSets;

% compute absolute cross-validation error using Statistics and Machine Learning Toolbox
% predfun = @(Atrain,ytrain,Atest)(Atest*fun(ytrain,Atrain));
% parfor i=1:n % for each component
%     err(i) = crossval('mse',A,y(:,i),'predfun',predfun,'partition',cvp);
% end

% compute relative cross-validation error
% if any(var(y,0,1)>eps)
%     ind = find(var(y,0,1)>eps);
%     err(ind) = err(ind)./var(y(:,ind),0,1); % err is divided by unbiaised empirical variance (second central moment) of y
% end
err = err./(moment(y,2,1)+mean(y,1).^2); % err is divided by empirical second moment of y

% compute corrected relative cross-validation error to reduce the sensitivity of error estimate to overfitting
% (non-corrected cross-validation error estimate underpredicts the error in L^2-norm (generalization error))
if ischarin('correction',varargin)
    C = inv(A'*A);
    if N~=P
        % corr = (N-1)/(N-P-1); % = (1-1/N)*(1-P/N-1/N)^(-1) % Adjusted Empirical Error (AEE) -> only accurate when N >> P (corr ~ 1+P/N when N->Inf)
        % corr = (N+P)/(N-P); % = (1+P/N)*(1-P/N)^(-1) % Future Prediction Error (FPE) [Akaike, 1970] -> only accurate when N >> P (corr ~ 1+2*P/N when N->Inf)
        % corr = max(0,(1-1*sqrt((P*(log(N/P)+1)-3)/N)))^(-1); % Uniform Convergence Bounds (UCB) [Cherkassky, Mulier, & Vapnik, 1997]
        corr = (N/(N-P))*(1+trace(C)); % = (1-P/N)^(-1)*(1+trace(C)) % Direct Eigenvalue Estimator [Chapelle, Vapnik & Bengio, 2002], [Blatman & Sudret, 2011] -> accurate even when N is not >> P (corr ~ 1+2*P/N when N->Inf, as trace(C) ~ P/N when N->Inf)
        err = err*corr;
    end
end

err = sqrt(err);

end
