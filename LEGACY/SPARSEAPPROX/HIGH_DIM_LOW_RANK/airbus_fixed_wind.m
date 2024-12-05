clear all
close all
warning off
fid = fopen('UC4_CHORUS-Input_fxWind.txt');
nCols = 28;
format = repmat('%f', [1 nCols]);
B = fscanf(fid, format,[nCols inf]);
fclose(fid);
A=B';
%%
fid = fopen('UC4_CHORUS-Output_fxWind.txt');
nCols = 6;
format = repmat('%f', [1 nCols]);
B = fscanf(fid, format,[nCols inf]);
fclose(fid);
FS=FS';
fsD=FS(:,1);
%%
xi = RANDVARS(RVUNIFORM(-1,1),28);
RV=RANDVARS(RVUNIFORM(0,1),28);
rsD(1:7000,1:28)=A(1:7000,1:28);
for i=1:28
    rsD(:,i) = transfer(RV{i},xi{i},rsD(:,i));
end
%%
fs=fsD(1:2000); fsModel=fsD(2001:3000);
rs=rsD(1:2000,:); rsModel=rsD(2001:3000,:);
%% Sparse least square
NIPARAM=struct('basis','global','order',1,'basistype','total');
[u,error]=sparse_least_square_solution(NIPARAM,fs,rs,fsModel,rsModel,RV);
%% Greedy
NIPARAM=struct('basis','global','basistype','total','order',2,'maxrank',10,'tol',...
    1.0e-03,'maxiter',20,'alphaupdate',0,'Kfold',1);
[u,error,dummy,dummy, dummy]=greedy_rank_one_solution(NIPARAM,fs,rs,fsModel,rsModel,xi);
%% Rank-M direct 
NIPARAM=struct('basis','global','basistype','total','order',2,'maxrank',10,'tol',...
    1.0e-03,'maxiter',20,'alphaupdate',0,'Kfold',1);
[u,error]=rank_m_solution(NIPARAM,fs,rs,fsModel,rsModel,xi);
%% Tensor Train using svd tolerance
NIPARAM=struct('basis','global','order',2,'svdtol',1.0e-08,'tol',...
            1.0e-06,'maxiter',15,'reg','l1');
[u,tt_error,rank,sr,err_norm]=tensor_train_solution_fast(NIPARAM,fs,rs,fsModel,rsModel,xi);
%% Tensor Train using rank selection with cross validation
NIPARAM=struct('basis','global','order',2,'svdmaxrank',10,'tol',...
            1.0e-06,'maxiter',15,'reg','l1');
[u,tt_error,rank,sr,err_norm]=tensor_train_solution_fast(NIPARAM,fs,rs,fsModel,rsModel,xi);
%% Tensor Train using rank adaptation
NIPARAM=struct('basis','global','order',2,'maxrank',10,'tol',...
            1.0e-06,'maxiter',15,'reg','l1');
[u,tt_error,rank,sr,err_norm]=tensor_train_solution_fast(NIPARAM,fs,rs,fsModel,rsModel,xi);