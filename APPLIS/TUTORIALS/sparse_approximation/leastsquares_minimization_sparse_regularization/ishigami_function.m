%% Least-squares minimization with sparse regularization - Ishigami function %%
%%---------------------------------------------------------------------------%%
% [Ishigami, Homma, 1990, ISUMA]
% [Saltelli, Chan, Scott, 2000, Wiley]
% [Sudret, 2008, RESS]
% [Blatman, Sudret, 2010, PEM]
% [Blatman, Sudret, 2010, RESS]
% [Blatman, Sudret, 2011, JCP]

% clc
clear all
close all

%% Random variables
M = 3; % number of random variables
rv = RVUNIFORM(-pi,pi);
RV = RANDVARS(repmat({rv},1,M));

%% Ishigami function
% y = sin(x_1) + a*sin(x_2)^2 + b*(x_3)^4*sin(x_1)
a = 7;
b = 0.1;
fun = @(x) (sin(x(:,1)) + a.*sin(x(:,2)).^2 + b.*x(:,3).^4.*sin(x(:,1)))';

%% Polynomial chaos basis
p = 9; % order of PC expansion
PC = POLYCHAOS(RV,p,'typebase',1); % PC basis
P = length(PC); % number of unknown PC expansion coefficients

%% Sampling
N = 300; % number of samples (regression points)
x = random(RV,N,1); % 1-by-M cell (associated to each random variable RV{1},...,RV{M})
                    % containing N-by-1 double (samples)
x = [x{:}]; % random samples (N-by-M matrix)
xpc = transfer(RV,RANDVARS(PC),x); % mapping of each random variable RV (defined over [-pi,pi])
                                   % into RANDVARS(PC) (defined over [-1,1] for Legendre polynomials)
                                   % for each random sample

%% Regularization
regul = 'l1'; % type of regularization ('' or 'l0' or 'l1')

%% Cross-validation
cv = 'leaveout'; % type of cross-validation procedure ('leaveout' or 'kfold')
k = 10; % number of folds

%% Least-squares minimization
t = tic;
y = fun(x)'; % evaluation of function fun on samples (N-by-1 vector)
A = polyval(PC,xpc); % evaluation of each polynomial of PC basis on samples (N-by-P matrix)
[u,err] = solve_leastsquares(A,y,regul,[],cv,k); % resolution of least-squares minimization problem
u = PCMATRIX(u,[1,1],PC); % solution (1-by-P vector)
time = toc(t);

%% Results
disp(' ')
disp(['M = ' num2str(getM(PC)) ' random variables'])
disp(['p = ' num2str(getorder(PC)) ' (order of PC expansion)'])
disp(['P = ' num2str(length(PC)) ' unknown PC expansion coefficients'])
disp(['N = ' num2str(N) ' samples (regression points)'])
disp(['N_opt = (M-1)*P = ' num2str((getM(PC)-1)*length(PC)) ' samples (regression points) according to empirical rule without regularization'])
disp(['eta = ' num2str(getSparsityRatio(u)) ' (sparsity index or ratio)'])
fprintf('error = %.4e (cross-validation error)\n',err)
fprintf('elapsed time = %f s\n',time);
disp(' ')

%% Quantities of interest : mean, variance, Sobol indices
% Analytical exact values
anal.mean = a/2;
anal.var = 1/2 + a^2/8 + b*pi^4/5 + b^2*pi^8/18;
anal.S1 = (1/2*(1 + b*pi^4/5)^2)/anal.var;
anal.S2 = (a^2/8)/anal.var;
anal.S3 = 0;
anal.S12 = 0;
anal.S13 = (b^2*pi^8*(8/225))/anal.var;
anal.S23 = 0;
anal.S123 = 0;
anal.S1T = anal.S1 + anal.S12 + anal.S13 + anal.S123;
anal.S2T = anal.S2 + anal.S12 + anal.S23 + anal.S123;
anal.S3T = anal.S3 + anal.S13 + anal.S23 + anal.S123;
% Numerical approximate values
num.mean = mean(u);
num.var = variance(u);
num.S1 = sobol_indices(u,1);
num.S2 = sobol_indices(u,2);
num.S3 = sobol_indices(u,3);
num.S12 = sobol_indices_group(u,[1,2]) - num.S1 - num.S2;
num.S13 = sobol_indices_group(u,[1,3]) - num.S1 - num.S3;
num.S23 = sobol_indices_group(u,[2,3]) - num.S2 - num.S3;
num.S123 = sobol_indices_group(u,[1,2,3]) - num.S1 - num.S2 - num.S3 - num.S12 - num.S13 - num.S23;
num.S1T = num.S1 + num.S12 + num.S13 + num.S123;
num.S2T = num.S2 + num.S12 + num.S23 + num.S123;
num.S3T = num.S3 + num.S13 + num.S23 + num.S123;
% Comparative table
fanal = '%10.5f';
fnum = '%9.5f';
ferr = '%14.4e';
disp('+-------------------+------------+-----------+----------------+')
disp('| Quantity \ Value  | Analytical | Numerical | Relative error |')
disp('+-------------------+------------+-----------+----------------+')
fprintf(['| Mean E            | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.mean,num.mean,abs((anal.mean-num.mean)/anal.mean))
fprintf(['| Variance V        | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.var,num.var,abs((anal.var-num.var)/anal.var))
fprintf(['| Sobol index S_1   | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S1,num.S1,abs((anal.S1-num.S1)/anal.S1))
fprintf(['| Sobol index S_2   | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S2,num.S2,abs((anal.S2-num.S2)/anal.S2))
fprintf(['| Sobol index S_3   | ' fanal ' | ' fnum ' |                |\n'],anal.S3,num.S3)
fprintf(['| Sobol index S_12  | ' fanal ' | ' fnum ' |                |\n'],anal.S12,num.S12)
fprintf(['| Sobol index S_13  | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S13,num.S13,abs((anal.S13-num.S13)/anal.S13))
fprintf(['| Sobol index S_23  | ' fanal ' | ' fnum ' |                |\n'],anal.S23,num.S23)
fprintf(['| Sobol index S_123 | ' fanal ' | ' fnum ' |                |\n'],anal.S123,num.S123)
fprintf(['| Sobol index S_1^T | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S1T,num.S1T,abs((anal.S1T-num.S1T)/anal.S1T))
fprintf(['| Sobol index S_2^T | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S2T,num.S2T,abs((anal.S2T-num.S2T)/anal.S2T))
fprintf(['| Sobol index S_3^T | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S3T,num.S3T,abs((anal.S3T-num.S3T)/anal.S3T))
disp('+-------------------+------------+-----------+----------------+')
