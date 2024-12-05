function [u,rank,sparsity_ratio,err_norm,rank_final,PC]=tensor_train_solution_fast(sparam,fs,rs,PC)

Hs=polyval(PC,rs);
Hs=cellfun(@full,Hs,'UniformOutput',false); 
dim=getnbgroups(PC);
p=getn(PC,dim);
n=size(fs,1);
% v=cell(1,dim);
rr=2;
r=rr*ones(1,dim+1);
r(1)=1;r(dim+1)=1;
if (strcmp(sparam.init,'canonical'))   
    [v,Tcan]=init_canonique(PC,fs,rs,p);
else
    for i=1:dim
        v{i}=rand(r(i),r(i+1),p);
    end
end

maxiter=sparam.maxiter
err_norm=zeros(maxiter,1);
rank=zeros(maxiter,dim-1);
sparsity_ratio=zeros(maxiter,dim);
for i=1:maxiter
    u=v;
    v = orth_tt(v);
    [E]=get_operator_matrices(v,dim,Hs);
    for k=1:dim-1
        
        alpha=getalpha2(E,k);
        beta=getbeta2(E,k+1);
        eta=geteta(alpha,beta,Hs,k);
        if isfield(sparam,'svdtol')
            A=permute(eta,[5 1 2 3 4]);
            s=size(eta);
            A=reshape(A,n,s(1)*s(2)*s(3)*s(4));
            if length(sparam.svdtol)~=1
                tolopt=get_opt_tolerance(A,fs,s,sparam.svdtol);
            else
                tolopt=sparam.svdtol;
            end
            param.lambda=0.;
            if (strcmp(sparam.reg,'ols')==1)
                 solopt=(A'*A)\(A'*fs);
            else
            [sol,path]=mexLasso(full(fs),full(A),param);
            %[solopt,epsil]=select_opt_path(A,fs,reduce_path(path,2),'leaveout','parallel');
            [solopt,epsil]=sel_opt_path_stat_test(A,fs,path);
            %sparsity_ratio(i,k)=nnz(solopt)/length(solopt);
            end
            W=reshape(solopt,s(1)*s(2),s(3)*s(4));
            [U,S,V]=svd(full(W));
            stt = diag(S);
            err=sqrt(1-cumsum(stt.^2)/sum(stt.^2));
            ko = find(err<tolopt);
            if isempty(ko)
                ko=numel(stt);
            else
                ko=min(ko);
            end
            
            %if ko>15
            %    ko=15;
            %end
            uu=reshape((U(:,1:ko)),size(v{k},1),p,ko);
            vv=reshape(V(:,1:ko)*S(1:ko,1:ko),size(v{k+1},2),p,ko);
            v{k}=permute(uu,[1 3 2]);
            v{k+1}=permute(vv,[3,1,2]);
        elseif isfield(sparam,'svdmaxrank')
            A=permute(eta,[5 1 2 3 4]);
            s=size(eta);
            A=reshape(A,n,s(1)*s(2)*s(3)*s(4));
            ko=get_opt_svd_rank(A,fs,s,sparam.svdmaxrank,sparam.reg);
            param.lambda=0.0;
            if (strcmp(sparam.reg,'ols')==1)
                 solopt=(A'*A)\(A'*fs);
            else
            [solopt,path]=mexLasso(full(fs),full(A),param);
            [solopt,epsil]=sel_opt_path_stat_test(A,fs,path);
            %[solopt,epsil]=select_opt_path(A,fs,path,'leaveout','correction');
            end
            sparsity_ratio(i,k)=nnz(solopt)/length(solopt);
            W=reshape(solopt,s(1)*s(2),s(3)*s(4));
            
            [U,S,V]=svd(full(W));
            if ko>numel(diag(S))
                ko=numel(diag(S));
            end
            uu=reshape((U(:,1:ko)),size(v{k},1),p,ko);
            vv=reshape(V(:,1:ko)*S(1:ko,1:ko),size(v{k+1},2),p,ko);
            v{k}=permute(uu,[1 3 2]);
            v{k+1}=permute(vv,[3,1,2]);
        else
            s=size(eta);        
            A=reshape(eta,s(1)*s(2),s(3)*s(4),s(5));
            [U,V,ko]=rankmin(A,fs,sparam.maxrank,sparam.reg);
 
uu=reshape(full(U),size(v{k},1),p,ko);
            vv=reshape(full(V),ko,size(v{k+1},2),p);
            v{k}=permute(uu,[1 3 2]);
            v{k+1}=vv;
            sparsity_ratio(i,:)=get_sparsity_ratio(v);
        end
        [E]=update_operator_matrices(E,v,k,Hs);
        rank(i,k)=ko;
 %       rank(1:i,:)
    end         
            rank_final=rank(i,:);
            e=tt_add(u,tt_minus(v));
            err_norm(i)=tt_norm(e)/tt_norm(u);
            fprintf('Iteration %d - Error %d\n',i,err_norm(i))
            fprintf('Rank = [%s]\n',num2str(rank_final))
            if err_norm(i)<sparam.tol %|| err_norm(i)>1.0e+03
                break
            end
end
%-------------------------------------------%
function [E]=get_operator_matrices(v,d,Hs)
E=cell(d,1);
    for j=1:d
        E{j}=ttm(v{j},Hs{j},3);
    end
%-------------------------------------------%
function [E]=update_operator_matrices(E,v,k,Hs)

    for j=k:k+1
        E{j}=ttm(v{j},Hs{j},3);
    end





%------------------------------------------------%
function [tolopt]=get_opt_tolerance(A,fs,s,svdtol)

[fsReg,AReg,fsTest,ATest]=cross_validation_Kfold_least_squares(A,fs,3);
nmodel=length(svdtol);
error=zeros(nmodel,3);
for m=1:nmodel
    tol=svdtol(m);
    for k=1:3
        param.lambda=0.;
        [sol,path]=mexLasso(full(fsReg{k}),full(AReg{k}),param);
        [solopt,epsil]=select_opt_path(AReg{k},fsReg{k},path,'leaveout');
        W=reshape(solopt,s(1)*s(2),s(3)*s(4));
        [U,Si,V]=svd(full(W));
        stt = diag(Si);
        err=sqrt(1-cumsum(stt.^2)/sum(stt.^2));
        tt = find(err<tol);
        if isempty(tt)
            tt=numel(stt);
        else
            tt=min(tt);
        end
        S=U(:,1:tt)*Si(1:tt,1:tt)*V(:,1:tt)';
        fapproxpath=ATest{k}*reshape(S,size(S,1)*size(S,2),1);
        error(m,k)=norm((fsTest{k}-fapproxpath),2)/norm(fsTest{k},2);
    end
    [eps c]=min(sum(error,2));
    tolopt=svdtol(c);
end
 
%-------------------------------------------------------------%
function [Uo,Vo,ko]=rankmin(A,fs,maxrank,reg)
M=3;
tt = min(maxrank,min(size(A,1),size(A ,2)));
[fsReg,AReg,fsTest,ATest]=cross_validation_Kfold_least_squares(A,fs,M);
error=zeros(tt,M);
for m=1:M

Ut=cell(1,tt);
Vt=cell(1,tt);
for k=1:tt
    Ut{k}=rand(size(A,1),k);
    [Ut{k},~]=qr(Ut{k},0);
    maxiter=3;
    for i=1:maxiter
        Vt{k}=minV_dmrg(AReg{m},Ut{k},fsReg{m},reg);
        [Vt{k},~] = qr(Vt{k}',0);
        Vt{k}=Vt{k}';
        Ut{k} =minU_dmrg(AReg{m},Vt{k},fsReg{m},reg);
        [Ut{k},R] = qr(Ut{k},0);
        Vt{k} = R*Vt{k};
        Snew = Ut{k}*Vt{k};
        if i>1
           errstagn = norm(S-Snew)/norm(S);
           if errstagn<1e-5
            break
           end
        end
        S = Snew;
    end
    S=Ut{k}*Vt{k};

    %sprintf('non zero factor = %f',nnz(S)/(size(S,1)*size(S,2)))
    BTest=permute(ATest{m},[3 1 2]);
    fapproxpath=reshape(BTest,size(BTest,1),size(BTest,2)*size(BTest,3))...
                *reshape(S,size(S,1)*size(S,2),1);
    error(k,m)=norm((fsTest{m}-fapproxpath),2)/norm(fsTest{m},2);
end
end

[~, ko]=min(sum(error,2));

Uo=Ut{ko};
Vo=Vt{ko};
%-----------------------------------------------------------
function [V] = minV_dmrg(A,U,b,reg)
B=ttm(A,U',1);
B=permute(B,[3 1 2]);
B=reshape(B,size(B,1),size(B,2)*size(B,3));
if (strcmp(reg,'ols')==1)
    Vo=(B'*B)\(B'*b);
else
param.lambda=0.0;
param.numThreads=-1;
[Vo,path]=mexLasso(full(b),full(B),param);
%[Vo]=select_opt_path(B,b,path,'leaveout','parallel');
[Vo]=sel_opt_path_stat_test(B,b,path);
end
if (nnz(Vo)==0)
    error('myApp:allzeros', 'All zeros selected in minV_dmrg')
end
V=reshape(full(Vo),size(U,2),size(A,2));
%-----------------------------------------------------------
function [U] = minU_dmrg(A,V,b,reg)
B=ttm(A,V,2);
B=permute(B,[3 1 2]);
B=reshape(B,size(B,1),size(B,2)*size(B,3));
if (strcmp(reg,'ols')==1)
    Uo=(B'*B)\(B'*b);
else
param.lambda=0.0;
param.numThreads=-1;
[Uo,path]=mexLasso(full(b),full(B),param);
%[Uo]=select_opt_path(B,b,path,'leaveout','parallel');
[Uo]=sel_opt_path_stat_test(B,b,path);
end
if (isempty(Uo))
    error('myApp:allzeros', 'All zeros selected in minU_dmrg')
end
U=reshape(full(Uo),size(A,1),size(V,1));
%-------------------------------------------------------------
function [U,V,tt]=initialise_rankmin(A,fs)
D=zeros(size(A,1),size(A,2));
for i =1:size(A,3)
    D=D+(A(:,:,i)*fs(i));
end
[U,Si,V]=svd(D);
%m=min(14,nnz(diag(Sigma)));
stt = diag(Si);
err=sqrt(1-cumsum(stt.^2)/sum(stt.^2));
tt = find(err<1e-16);
if isempty(tt)
    tt=numel(stt);
else
    tt=min(tt);
end
%--------------------------------------------------------
function [w]=tt_add(u,v)
    d=size(u,2);
    w=cell(1,size(u,2));
    w{1}=horzcat(u{1},v{1});
    w{d}=vertcat(u{d},v{d});
    for i=2:d-1
        a=cell(2);
        a{1,1}=u{i};
        a{1,2}=zeros(size(u{i},1),size(v{i},2),size(u{i},3));
        a{2,1}=zeros(size(v{i},1),size(u{i},2),size(u{i},3));
        a{2,2}=v{i};
        w{i}=cell2mat(a);
    end
%-----------------------------------------------------------
function [w]=tt_minus(v)
v{1}=-1*v{1};
w=v;
%----------------------------------------------------------
function [w]=tt_dot(A,B)
w=1;
d=size(A,2);
for i=1:d
    W{i}=cell(size(A{i},1),size(A{i},2));
for j=1:size(A{i},1)
    for k=1:size(A{i},2)
        u=zeros(size(A{i},3),1);
        u(:)=A{i}(j,k,:);
        W{i}{j,k}=ttv(B{i},u,3);
    end
end
W{i}=cell2mat(W{i});
end
for i=1:d
    w=w*W{i};
end
%------------------------------------------------------------
function [w]=tt_norm(A)
w=sqrt(tt_dot(A,A));
%----------------------------------------------------------
% function [alpha]=getalpha1(v,k,Hs,n)
% alpha=cell(n,1);
% alpha(:)={1};
% A=cell(1,k-1);
% for ii=1:n
% for j=1:k-1
% A{j}=ttv(v{j},Hs{j}(ii,:)',3);
% end
% for i=1:k-1
%     alpha{ii}=alpha{ii}*A{i};
% end
% end
%----------------------------------------------------------
function [alpha]=getalpha2(E,k)
%alpha=cell(n,1);
%alpha(:)={1};
n = size(E{1},3);
alpha=ones(n,1);
for mu=1:k-1
    alpha = repmat(alpha,[1 1 size(E{mu},2)]);
    alpha = sum(alpha.*permute(E{mu},[3 1 2]),2);
    alpha = alpha(:,:);
end

return
%----------------------------------------------------------
function [beta]=getbeta2(E,k)
d=numel(E);
n = size(E{1},3);
beta=ones(n,1);
for mu=d:-1:k+1
    beta = repmat(beta,[1 1 size(E{mu},1)]);
    beta = sum(beta.*permute(E{mu},[3 2 1]),2);
    beta = beta(:,:);
end

return

%------------------------------------------------------------
% function [beta]=getbeta1(v,k,Hs,n)
% beta=cell(n,1);
% beta(:)={1};
% d=size(v,2);
% A=cell(1,d-k);
% for ii=1:n
%     for j=k+1:d
%         A{j}=ttv(v{j},Hs{j}(ii,:)',3);
%     end
%     for i=k+1:d
%         beta{ii}=beta{ii}*A{i};
%     end
% end
%------------------------------------------------------------
function [eta]=geteta(alpha,beta,Hs,k)
n = size(Hs{k},1);
s = [size(alpha,2),size(Hs{k},2),size(beta,2),size(Hs{k+1},2),n];

alpha = permute(alpha,[2 3 4 5 1]);
alpha = repmat(alpha,[1 s(2) s(3) s(4) 1]);
beta = permute(beta,[3 4 2 5 1]);
beta = repmat(beta,[s(1) s(2) 1 s(4) 1]);
H1 = permute(Hs{k},[3  2 4 5 1]);
H1 = repmat(H1,[s(1) 1 s(3) s(4) 1]);
H2 = permute(Hs{k+1},[3 4 5 2 1]);
H2 = repmat(H2,[s(1) s(2) s(3) 1 1]);

eta = alpha.*beta.*H1.*H2;

return
%------------------------------------------------------------
function [s]=get_sparsity_ratio(v)
d=length(v);
s=zeros(1,d);
for i=1:d
    s(i)=nnz(v{i})/numel(v{i});
end
%------------------------------------------------------------
function [redpath]=reduce_path(path,d)
npath=size(path,2);
rpath=1:d:npath;
redpath=zeros(size(path,1),size(rpath,2));
redpath=path(:,rpath);

