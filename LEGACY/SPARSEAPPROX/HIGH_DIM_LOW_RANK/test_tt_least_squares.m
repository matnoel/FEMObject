%clear all
cas=3;
switch cas
    case {2,5}
        fun = @(x) multivariateFunctionsSamplesTestcases(cas,x);  
case 1
    fun = @(x) multivariateFunctionsSamplesTestcases(1,x,5);
case 3
    d=5;w = rand();c=rand(1,d);
fun = @(x) multivariateFunctionsSamplesTestcases(3,x,d,w,c);  
case 4
    d=8;w = rand(1,d);c=rand(1,d);
fun = @(x) multivariateFunctionsSamplesTestcases(4,x,d,w,c);  

end
%%
close all
[rsModel,fsModel,RV]=fun(1000);

%%
ssize=[100];
deg=[3];
ls_tot_error=zeros(length(ssize),length(deg));
ls_partial_error=zeros(length(ssize),length(deg));
greedy_error=zeros(length(ssize),length(deg));
rank_m_error=zeros(length(ssize),length(deg));
tt_error=zeros(length(ssize),length(deg));
for s=1:length(ssize)
[rs,fs,RV]=fun(ssize(s));
for pp=1:length(deg)
    clear NIPARAM
NIPARAM=struct('basis','global','order',deg(pp),'init','canonical',...
            'tol',1.0e-08,'maxiter',5,'alphaupdate',0,...
            'reg','l1');
%NIPARAM =setfield(NIPARAM,'svdmaxrank',10);
NIPARAM =setfield(NIPARAM,'maxrank',8);
%NIPARAM =setfield(NIPARAM,'svdtol',1e-2);
[X,PC]=PCTPMODEL(RV,'order',deg(pp),'pcg','typebase',1,'nomasse');
[u,rank,sr]=tensor_train_solution_fast(NIPARAM,fs,rs,PC);
HsTest=polyval(PC,rsModel);
fsTest=evaluateTT(u,HsTest);
tt_error(s,pp)=norm(fsModel-fsTest,2)/norm(fsModel,2);
end
end
tt_error

%%
[rsModel,fsModel,RV]=fun(1000);

HsTest=polyval(PC,rsModel);
fsTest=evaluateTT(u,HsTest);
modelerror=norm(fsModel-fsTest,2)/norm(fsModel,2)


