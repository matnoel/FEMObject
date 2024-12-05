function [u,modelerror,optrank,delta_tot,rank_modelerror]=greedy_rank_one_solution(param,fs,rs,fsModel,rsModel,RV)

if (strcmp(param(1).basis,'global'))   
    [X,PC]=PCTPMODEL(RV,'order',param.order,'pcg','typebase',1,'nomasse');
elseif (strcmp(param(1).basis,'fe'))
    [X PC]=PCTPMODEL(RV,'order',param.order,'pcg','fedim',1:2,'femesh',{8,8},'typebase',1,'nomasse');
elseif (strcmp(param(1).basis,'wavelet'))
    [X PC] = PCTPMODEL(RV,'order',param.order,'pcg','waveletdim',1:2,'typebase',1,'waveletlevel',param.waveletlevel,'nomasse');
end
if param.Kfold>1
[fsReg,rsReg,fsTest,rsTest]=cross_validation(fs,rs,param.Kfold);  
error=zeros(param(1).maxrank,1);
for k=1:param.Kfold
    utemp=greedy_construction(fsReg{k},rsReg{k},PC,param(1).maxrank,param(1).tol,param(1).alphaupdate,param(1).maxiter);
    errortemp=get_greedy_construction_error(utemp,fsTest{k},rsTest{k},PC,param(1).maxrank);
    error=error+errortemp;
end
[c rankopt] = min(abs(error));
[u,delta_tot]=greedy_construction(fs,rs,PC,rankopt,param(1).tol,param(1).alphaupdate,param(1).maxiter);
else
[u,delta_tot]=greedy_construction(fs,rs,PC,param.maxrank,param(1).tol,param(1).alphaupdate,param(1).maxiter);    
end
rank_modelerror=get_greedy_construction_error(u,fsModel,rsModel,PC,getm(u));
[modelerror optrank]=min(rank_modelerror);
%---------------------------------%
function [utemp,delta_tot]=greedy_construction(fstemp,rstemp,PC,maxrank,tol,alphaupdate,maxiter)
utemp=PCTPMATRIXSUM(PC);
fsm=fstemp;
delta=zeros(length(fstemp),maxrank);
for m=1:maxrank
[lopt,delta(:,m)]=get_optimal_l(PC,fsm,rstemp,tol,maxiter);
us=evaluate_rank_one_fn(PC,lopt,rstemp);
if isempty(us)
   fprintf('Max rank is %d\n',m-1);
   break
else
utemp=utemp+lopt;
if (alphaupdate)
    [u_updated coeff]=alpha_update(fstemp,rstemp,PC,utemp);
    utemp=u_updated;
end
fsm=calc_residual(PC,utemp,fstemp,rstemp);
end
end
delta_tot=sum(delta,2);
%--------------------------------%
function [errortemp]=get_greedy_construction_error(u,fsTest,rsTest,PC,maxrank)
errortemp=zeros(maxrank,1);
modelapprox=zeros(size(fsTest,1),1);
for m=1:getm(u)
    modelapprox=modelapprox+evaluate_rank_one_fn(PC,u{m},rsTest);
    errortemp(m)=norm((fsTest-modelapprox),2)/norm(fsTest,2);
end
%--------------------------------%
function [fsReg,rsReg,fsTest,rsTest]=cross_validation(fs,rs,K)
fsReg={};rsReg={}; fsTest={}; rsTest={};
CVO=cvpartition(length(fs),'kfold',K);
for i=1:K
    trIdx=CVO.training(i);
    fsReg{i}=fs(trIdx); rsReg{i}=rs(trIdx,:);
    fsTest{i}=fs(~trIdx); rsTest{i}=rs(~trIdx,:);
end
%-----------------------------%
function [lopt,delta]=get_optimal_l(PC,fstemp,rstemp,tol,maxiter)
delta=zeros(length(fstemp),1);
dim=getnbgroups(PC);
Hs=polyval(PC,rstemp);
l=normalize(normalizephi(ones(1,1,PC)));
ls=cell(1,dim);
for i=1:dim
    ls{i} = Hs{i}*l{i};
end
for pp=1:maxiter
    lp=l;
    for i=1:dim       
        l{0}=1;  
        ls{i}=ones(length(fstemp),1);
        gammas = prod([ls{:}],2);
        D = bsxfun(@times, Hs{i}, gammas);
%         Dred=D(:,any(D)); 
%         nzbasis= any(D)~=0;
        param.lambda=1.0e-03;
        param.numThreads=-1; 
        path=[];     
        [sol path]=mexLasso(full(fstemp),full(D),param);
        %[solred path]=mexLasso(full(fstemp),full(Dred),param);
        if nnz(path)>0
            [sol,dummy1,delta]=sel_opt_path_stat_test(D,fstemp,path);
            %[sol,dummy1,dummy2,delta]=select_opt_path(D,fstemp,path,'leaveout','correction');
%             sol=zeros(size(D,2),1);
%             if nnz(solred)
%                 sol(nzbasis)=solred;
%             else
%                 break;
%             end 
        else
            break;
        end
        l{i}=sol;
        normphi = norm(l{i});
        l{i} = l{i}/normphi;
        ls{i}= Hs{i}*l{i};
        l{0}=normphi;
    end
    error=norm(lp-l);
    fprintf('iter %d - err %d\n',pp,error);
    if(error<tol)
        break;
    end
end
lopt=l;
%-----------------------------%
% function [pathOpt epsil]=sel_opt_path(D,fstemp,path)
% path(:,1)=[];
% path = evaluate_ls_coeff(D,fstemp,path);
% epsil=zeros(size(path,2),1);
% for i = 1:size(path,2)
%     fapproxpath=D*path(:,i);
%     ind=find(path(:,i));
%     P=D(:,[ind]);
%     G=(P'*P)^(-1);
%     C=P*G*P';
%     merror=((fstemp-fapproxpath)./(1-diag(C))).^2;
%     error=sum(merror,1)/size(fstemp,1);
%     epsil(i) = error/var(fstemp);
%     corrFact=(1/(1-nnz(path(:,1))/size(fstemp,1)))*(1+trace(G));
%     epsil(i)=epsil(i)*corrFact;
% end
% [c j] = min(abs(epsil));
% pathOpt = path(:,j);
%-----------------------------------%
% function [path_ls]=evaluate_ls_coeff(D,fstemp,path)
% path_ls=zeros(size(path));
% for i=1:size(path,2)
%     sol=path(:,i);
%     ind=find(sol);
%     Dred=D(:,[ind]);
%     x=(Dred'*Dred)\(Dred'*fstemp);
%     vec=zeros(size(sol,1),1);
%     vec(ind)=x;
%     path_ls(:,i)=vec;
% end
%-------------------------%
function [u_updated coeff]=alpha_update(fs,rs,PC,u)
    l = getfuns(u);
    for jj=1:getm(u);
    l{jj}{0}=1;
    end
    for j=1:getm(u)
    tensorBasisAll(:,j)=evaluate_rank_one_fn(PC,l{j},rs);    
    end
    fsReg = {};tensorBasisReg={}; fsTest = {}; tensorBasisTest = {};
    K=min(3,length(fs));
    CVO=cvpartition(length(fs),'kfold',K);
    for i=1:K
    trIdx = CVO.training(i);
    fsReg{i} = fs(trIdx); tensorBasisReg{i}=tensorBasisAll(trIdx,:);
    fsTest{i} = fs(~trIdx); tensorBasisTest{i}=tensorBasisAll(~trIdx,:);
    end

    for k=1:K
    A=tensorBasisReg{k};
    b=fsReg{k};     
    testBasis=tensorBasisTest{k};
    param.lambda=0.0;        
    param.numThreads=-1; 
    [alpha path]=mexLasso(full(b),full(A),param);        
    [pathOpt(:,k)]=calc_cverror(fsTest{k},path,testBasis);
    end
    sol=calc_cverror(fs,pathOpt,tensorBasisAll);
    for jj=1:getm(u)
    l{jj}{0}=sol(jj);
    end
    u_updated = PCTPMATRIXSUM(l{:});
    coeff=sol;
%---------------------%
function [pathOpt ind]=calc_cverror(fsTest,path,testBasis)
epsil=zeros(size(path,2),1);
for i = 1:size(path,2)
    error = fsTest-(testBasis*path(:,i));
    epsil(i) = (norm(error,2)/norm(fsTest,2))/length(fsTest);
end
[c j] = min(epsil);
pathOpt = path(:,j);
ind=j;
%---------------------%
function [residual]=calc_residual(PC,u,fs,rs)
modelapprox=zeros(size(rs,1),1);
for m=1:getm(u)
    modelapprox=modelapprox+evaluate_rank_one_fn(PC,u{m},rs);
end
residual=fs-modelapprox;
%----------------------%
function [feval]=evaluate_rank_one_fn(PC,lopt,rs)
if (isempty(getphi0(lopt)))
        feval=[];
else
modelTestBasis=polyval(PC,rs);
modelapprox=ones(size(rs,1),1)*lopt{0};
for i=1:getnbgroups(PC)
    modelapprox=modelapprox.*(modelTestBasis{i}*lopt{i});
end
feval=modelapprox;
end
%-----------------------%



    