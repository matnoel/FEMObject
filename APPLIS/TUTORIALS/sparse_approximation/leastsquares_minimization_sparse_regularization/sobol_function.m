%% Least-squares minimization with sparse regularization - Sobol function %%
%%------------------------------------------------------------------------%%
% [Saltelli, Sobol, RESS, 1995]
% [Sobol, 2003, RESS]
% [Sudret, 2008, RESS]

% clc
clear all
close all

%% Random variables
M = 8; % number of random variables
rv = RVUNIFORM(0,1);
RV = RANDVARS(repmat({rv},1,M));

%% Sobol function - polynomial function of degree M = 8
% y = prod_{j=1}^{M}(|4*x_j-2|+a_j)/(1+a_j)
a = [1,2,5,10,20,50,100,500];
% fun = @(x) ((abs(4*x(:,1)-2)+a(1))/(1+a(1)).*(abs(4*x(:,2)-2)+a(2))/(1+a(2)).*...
%     (abs(4*x(:,3)-2)+a(3))/(1+a(3)).*(abs(4*x(:,4)-2)+a(4))/(1+a(4)).*...
%     (abs(4*x(:,5)-2)+a(5))/(1+a(5)).*(abs(4*x(:,6)-2)+a(6))/(1+a(6)).*...
%     (abs(4*x(:,7)-2)+a(7))/(1+a(7)).*(abs(4*x(:,8)-2)+a(8))/(1+a(8)))';
fun = @(x) (prod((abs(4*x-2)+repmat(a,size(x,1),1))./repmat(1+a,size(x,1),1),2))';

%% Polynomial chaos basis
p = 2; % order of PC expansion
PC = POLYCHAOS(RV,p,'typebase',1); % PC basis
P = length(PC); % number of unknown PC expansion coefficients

%% Sampling
N = 200; % number of samples (regression points)
x = random(RV,N,1); % 1-by-M cell (associated to each random variable RV{1},...,RV{M})
                    % containing N-by-1 double (samples)
x = [x{:}]; % random samples (N-by-M matrix)
xpc = transfer(RV,RANDVARS(PC),x); % mapping of each random variable RV (defined over [0,1])
                                   % into RANDVARS(PC) (defined over [-1,1] for Legendre polynomials)
                                   % for each random sample

%% Regularization
regul = 'l1'; % type of regularization ('' or 'l0' or 'l1')

%% Cross-validation
cv = 'leaveout'; % type of cross-validation procedure ('leaveout' or 'kfold')
k = 10; % number of folds (only for k-fold cross-validation procedure)

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
anal.mean = 1;
anal.var = 1/(1+a(1))^2*(4/3+2*a(1)+a(1)^2)*1/(1+a(2))^2*(4/3+2*a(2)+a(2)^2)*1/(1+a(3))^2*(4/3+2*a(3)+a(3)^2)...
    *1/(1+a(4))^2*(4/3+2*a(4)+a(4)^2)*1/(1+a(5))^2*(4/3+2*a(5)+a(5)^2)*1/(1+a(6))^2*(4/3+2*a(6)+a(6)^2)...
    *1/(1+a(7))^2*(4/3+2*a(7)+a(7)^2)*1/(1+a(8))^2*(4/3+2*a(8)+a(8)^2) - 1;
anal.S1 = 1/anal.var*1/(3*(1+a(1))^2);
anal.S2 = 1/anal.var*1/(3*(1+a(2))^2);
anal.S3 = 1/anal.var*1/(3*(1+a(3))^2);
anal.S4 = 1/anal.var*1/(3*(1+a(4))^2);
anal.S5 = 1/anal.var*1/(3*(1+a(5))^2);
anal.S6 = 1/anal.var*1/(3*(1+a(6))^2);
anal.S7 = 1/anal.var*1/(3*(1+a(7))^2);
anal.S8 = 1/anal.var*1/(3*(1+a(8))^2);
anal.S12 = 1/anal.var*1/(3*(1+a(1))^2)*1/(3*(1+a(2))^2);
anal.S13 = 1/anal.var*1/(3*(1+a(1))^2)*1/(3*(1+a(3))^2);
anal.S14 = 1/anal.var*1/(3*(1+a(1))^2)*1/(3*(1+a(4))^2);
anal.S15 = 1/anal.var*1/(3*(1+a(1))^2)*1/(3*(1+a(5))^2);
anal.S16 = 1/anal.var*1/(3*(1+a(1))^2)*1/(3*(1+a(6))^2);
anal.S17 = 1/anal.var*1/(3*(1+a(1))^2)*1/(3*(1+a(7))^2);
anal.S18 = 1/anal.var*1/(3*(1+a(1))^2)*1/(3*(1+a(8))^2);
% Numerical approximate values
num.mean = mean(u);
num.var = variance(u);
num.S1 = sobol_indices(u,1);
num.S2 = sobol_indices(u,2);
num.S3 = sobol_indices(u,3);
num.S4 = sobol_indices(u,4);
num.S5 = sobol_indices(u,5);
num.S6 = sobol_indices(u,6);
num.S7 = sobol_indices(u,7);
num.S8 = sobol_indices(u,8);
num.S12 = sobol_indices_group(u,[1,2]) - num.S1 - num.S2;
num.S13 = sobol_indices_group(u,[1,3]) - num.S1 - num.S3;
num.S14 = sobol_indices_group(u,[1,4]) - num.S1 - num.S4;
num.S15 = sobol_indices_group(u,[1,5]) - num.S1 - num.S5;
num.S16 = sobol_indices_group(u,[1,6]) - num.S1 - num.S6;
num.S17 = sobol_indices_group(u,[1,7]) - num.S1 - num.S7;
num.S18 = sobol_indices_group(u,[1,8]) - num.S1 - num.S8;
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
fprintf(['| Sobol index S_3   | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S3,num.S3,abs((anal.S3-num.S3)/anal.S3))
fprintf(['| Sobol index S_4   | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S4,num.S4,abs((anal.S4-num.S4)/anal.S4))
fprintf(['| Sobol index S_5   | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S5,num.S5,abs((anal.S5-num.S5)/anal.S5))
fprintf(['| Sobol index S_6   | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S6,num.S6,abs((anal.S6-num.S6)/anal.S6))
fprintf(['| Sobol index S_7   | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S7,num.S7,abs((anal.S7-num.S7)/anal.S7))
fprintf(['| Sobol index S_8   | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S8,num.S8,abs((anal.S8-num.S8)/anal.S8))
fprintf(['| Sobol index S_12  | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S12,num.S12,abs((anal.S12-num.S12)/anal.S12))
fprintf(['| Sobol index S_13  | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S13,num.S13,abs((anal.S13-num.S13)/anal.S13))
fprintf(['| Sobol index S_14  | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S14,num.S14,abs((anal.S14-num.S14)/anal.S14))
fprintf(['| Sobol index S_15  | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S15,num.S15,abs((anal.S15-num.S15)/anal.S15))
fprintf(['| Sobol index S_16  | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S16,num.S16,abs((anal.S16-num.S16)/anal.S16))
fprintf(['| Sobol index S_17  | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S17,num.S17,abs((anal.S17-num.S17)/anal.S17))
fprintf(['| Sobol index S_18  | ' fanal ' | ' fnum ' | ' ferr ' |\n'],anal.S18,num.S18,abs((anal.S18-num.S18)/anal.S18))
disp('+-------------------+------------+-----------+----------------+')

% Results
% Only the three or four first variables have a significant influence on the solution variance
