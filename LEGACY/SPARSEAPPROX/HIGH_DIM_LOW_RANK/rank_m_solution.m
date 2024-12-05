function [u,modelerror,error]=rank_m_solution(param,fs,rs,fsModel,rsModel,RV)
if (strcmp(param.basis,'global'))
    [X,PC]=PCTPMODEL(RV,'order',param.order,'pcg','typebase',1,'nomasse');
elseif (strcmp(param.basis,'fe'))
    [X PC]=PCTPMODEL(RV,'order',param.order,'pcg','fedim',1:2,'femesh',{6,6},'typebase',1,'nomasse');
elseif (strcmp(param.basis,'wavelet'))
    [X PC]=PCTPMODEL(RV,'order',param.order,'pcg','waveletdim',1:2,'typebase',1,'waveletlevel',{3,3},'nomasse');
end
[fsReg,rsReg,fsTest,rsTest]=cross_validation(fs,rs,3);  
error=zeros(param.maxrank,1);
for rank=1:param.maxrank
    for k=1:param.Kfold
        utemp=rank_m_construction(fsReg{k},rsReg{k},PC,rank,param,fsTest{k},rsTest{k});
        error(rank)=error(rank)+get_rank_m_construction_error(utemp,fsTest{k},rsTest{k},PC,rank);
    end
end
[c rankopt] = min(abs(error));
u=rank_m_construction(fs,rs,PC,rankopt,param,fsModel,rsModel);
modelerror=get_rank_m_construction_error(u,fsModel,rsModel,PC,getm(u));
%----------------------------------%
function [u]=rank_m_construction(fs,rs,PC,rank,sparam,fsTest,rsTest)
m=rank;
p=getn(PC,1);
d=size(rs,2);
F=cell(m,d);
for i=1:m
    for j=1:d
        F{i,j}= speye(p,1);
    end
end
alpha=ones(m,1);
Hs=polyval(PC,rs);
A=cell(1,m);
if (strcmp(sparam.reg,'l2'))
    HsTest=polyval(PC,rsTest);
    ATest=cell(1,m);
end
l=SEPMATRIX(F,alpha);
for k=1:sparam.maxiter
    l0 = l;
for i=1:d
    for j=1:m
        gamma=getgamma(Hs,l.F,j,i,size(fs,1));
        if (strcmp(sparam.reg,'l2'))
            gammaTest=getgamma(HsTest,l.F,j,i,size(fsTest,1));
        end
        A{j}=(Hs{i}).*(repmat(gamma,1,size(Hs{i},2)));
        if (strcmp(sparam.reg,'l2'))
            ATest{j}=(HsTest{i}).*(repmat(gammaTest,1,size(HsTest{i},2)));    
        end
    end
    D=horzcat(A{:});
    if (strcmp(sparam.reg,'l2'))
        DTest=horzcat(ATest{:});
    end
    param.lambda=0;
    param.numThreads=-1;
    if (strcmp(sparam.reg,'l1'))
        [sol,path]=mexLasso(full(fs),full(D),param);
        [sol,epsil]=sel_opt_path_stat_test(D,fs,path);
    elseif (strcmp(sparam.reg,'l2'))
        sol=get_opt_l2_sol(D,fs,p,m,DTest,fsTest);
        lambda=1.0e-02;
        %sol=(D'*D+lambda*eye((p)*m))\(D'*fs);
    else
        sol=(D'*D)\(D'*fs);
    end
    if(isempty(sol))
        fprintf('All zero coefficients selected, Matrix coherence:%f \n',matrix_coherence(D));
        break
    end
    sol=reshape(sol,p,m);
    for j=1:m
        if all(sol(:,j)==0)
            l.F{j,i}=sol(:,j);
        else
            l.alpha(j)=norm(sol(:,j));
            l.F{j,i}=sol(:,j)/l.alpha(j);
        end
    end
end
err = norm(PCTPMATRIX(l0-l,PC))/norm(PCTPMATRIX(l0,PC)); % gives error in mulitplyF
fprintf('iter %d - err %d\n',k,err);
if err<sparam.tol 
    break
end
end
u=PCTPMATRIX(l,PC);
%-----------------------------------%
function [fsReg,rsReg,fsTest,rsTest]=cross_validation(fs,rs,K)
fsReg={};rsReg={}; fsTest={}; rsTest={};
CVO=cvpartition(length(fs),'kfold',K);
for i=1:K
    trIdx=CVO.training(i);
    fsReg{i}=fs(trIdx); rsReg{i}=rs(trIdx,:);
    fsTest{i}=fs(~trIdx); rsTest{i}=rs(~trIdx,:);
end
%------------------------------------%
function [gamma]=getgamma(Hs,F,j,i,numSample)
d=size(Hs,2);
x=zeros(numSample,d);
for k=1:d
    x(:,k)=(Hs{k}*F{j,k});
end
x(:,i)=ones(numSample,1);
gamma=prod(x,2);  
%-------------------------------------%
function [error]=get_rank_m_construction_error(u,fsTest,rsTest,PC,maxrank)
modelapprox=zeros(size(fsTest,1),1);
for m=1:getm(u)
    if(strcmp(class(u),'PCTPMATRIX'))
        modelapprox=modelapprox+evaluate_rank_one_fn(PC,u,rsTest);
    else
        modelapprox=modelapprox+evaluate_rank_one_fn(PC,u{m},rsTest);
    end
end
error=norm((fsTest-modelapprox),2)/norm(fsTest,2);
%--------------------------------------%
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
%----------------------------------------%
function [pathOpt epsil]=sel_opt_path(D,fstemp,path)
path(:,1)=[];
path = evaluate_ls_coeff(D,fstemp,path);
epsil=zeros(size(path,2),1);
for i = 1:size(path,2)
    fapproxpath=D*path(:,i);
    ind=find(path(:,i));
    P=D(:,[ind]);
    G=(P'*P)^(-1);
    C=P*G*P';
    merror=((fstemp-fapproxpath)./(1-diag(C))).^2;
    error=sum(merror,1)/size(fstemp,1);
    epsil(i) = error/var(fstemp);
    corrFact=(1/(1-nnz(path(:,1))/size(fstemp,1)))*(1+trace(G));
    epsil(i)=epsil(i)*corrFact;
end
[c j] = min(abs(epsil));
pathOpt = path(:,j);
%------------------------------------------%
function [pathOpt]=get_opt_l2_sol(D,fs,p,m,DTest,fsTest)
Lambda=[1.0e-01,1.0e-02,1.0e-03,1.0e-04,1.0e-05,1.0e-06,1.0e-07];
L=length(Lambda);
path=zeros(size(D,2),L);
error=zeros(L,1);
for i=1:length(Lambda)
    path(:,i)=(D'*D+Lambda(i)*eye((p)*m))\(D'*fs);
    fapproxpath=DTest*path(:,i);
    error(i)=norm((fsTest-fapproxpath),2);
end
[c j] = min(error);
pathOpt = path(:,j);
%-----------------------------------%
function [path_ls]=evaluate_ls_coeff(D,fstemp,path)
path_ls=zeros(size(path));
for i=1:size(path,2)
    sol=path(:,i);
    ind=find(sol);
    Dred=D(:,[ind]);
    x=(Dred'*Dred)\(Dred'*fstemp);
    vec=zeros(size(sol,1),1);
    vec(ind)=x;
    path_ls(:,i)=vec;
end
%------------------------------------%
function [coherence]=matrix_coherence(A)
B = zeros(size(A));
for i=1:size(A,2)
    for j=1:size(A,2)
        if j>i
            B(i,j)=abs(A(:,i)'*A(:,j))/(norm(A(:,i),2)*norm(A(:,j),2));
        else
            B(i,j)=0;
        end
    end
end
coherence=max(max(B));