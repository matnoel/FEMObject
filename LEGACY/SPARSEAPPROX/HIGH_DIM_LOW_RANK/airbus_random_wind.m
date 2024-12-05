clear all
close all
fid = fopen('UC4_CHORUS-Input_fxWind.txt');
B=fscanf(fid, '%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g', ...
                [29 inf]);
fclose(fid);
A=B';
%%
fid = fopen('UC4_CHORUS-Output_fxWind.txt');
FS=fscanf(fid, '%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g', ...
                [6 inf]);
fclose(fid);
FS=FS';
fsD=FS(:,1);
%%
xi = RANDVARS(RVUNIFORM(-1,1),29);
RV=RANDVARS(RVUNIFORM(0,1),29);
rsD(1:7000,1:29)=A(1:7000,1:29);
for i=1:29
    rsD(:,i) = transfer(RV{i},xi{i},rsD(:,i));
end
%%
fs=fsD(1:2000); fsModel=fsD(2001:5000);
rs=rsD(1:2000,:); rsModel=rsD(2001:5000,:);
%%
% deg=[2,3,4,5];
% modelerror=zeros(length(deg),1);
% modelerror_rank_m=zeros(length(deg),1);
% greedy_error=cell(length(deg),1);
% rank_m_error=cell(length(deg),1);
% greedy_u=cell(length(deg),1);
% rank_m_u=cell(length(deg),1);
% for pp=1:length(deg)
% NIPARAM=struct('basis','global','basistype','total','order',deg(pp),'maxrank',10,'tol',...
%     1.0e-08,'maxiter',20,'alphaupdate',0,'Kfold',1);
% [greedy_u{pp},modelerror(pp,1),dummy,dummy, greedy_error{pp}]=greedy_rank_one_solution(NIPARAM,fs,rs,fsModel,rsModel,xi);
% [rank_m_u{pp},modelerror_rank_m(pp,1),rank_m_error{pp}]=rank_m_solution(NIPARAM,fs,rs,fsModel,rsModel,xi);
% end
% %%
% S=diag(sobol_indices(u));
% for ii=1:29
%     for jj=1:ii-1
%         S(ii,jj)=calc_sobol(u,[ii,jj])-S(ii,ii)-S(jj,jj);
%     end
% end
% imagesc(S);   
% %%
degree=3;
K=[5,7];
lambda=[1.0e-06];
NIPARAM=struct('basis','global','order',degree,'maxrank',2,'tol',...
    1.0e-03,'maxiter',10,'alphaupdate',0,'Kfold',1);
[clust_error,lambda_error,err_sans_clust,final_cluster]=lr_classification(fs,rs,fsModel,rsModel,RV,NIPARAM,K,lambda);
% %%
% NIPARAM=struct('basis','global','order',4,'svdtol',1.0e-05,'tol',...
%             1.0e-06,'maxiter',20,'reg','l1');
% [u,tt_error,rank,sr]=tensor_train_solution_fast(NIPARAM,fs,rs,fsModel,rsModel,xi);