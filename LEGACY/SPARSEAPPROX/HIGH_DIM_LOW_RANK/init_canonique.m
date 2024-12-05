function [v,u]=init_canonique(PC,fs,rs,p)
u=PCTPMATRIXSUM(PC);
maxRank=2;
fsm=fs;
for m=1:maxRank
    fprintf('m = %d\n',m)
    [lopt]=get_optimal_l(PC,fsm,rs);
    u=u+lopt;
    us=evaluate_rank_one_fn(PC,lopt,rs);
    if isempty(us)
        fprintf('Max rank is %d\n',m-1);
        break
    else
        fsm=fsm-us;
    end
end
maxRank = m-1;
d=getnbdim(PC);
% v=cell(1,d);
% for j=1:5
% v{1}(1,j,:)=double(u{j}{1});
% v{d}(j,1,:)=double(u{j}{d});
% end
% for i=2:d-1
%     for j=1:maxrank
%     v{i}=vertcat(v{i},u{j}{i});
%     end
%     v{i}=reshape(v{i},[p+1,5,5]);
%     v{i}=permute(v{i},[2,3,1]);
% end
v = cell(1,d);
v{1} = zeros(1,maxRank,p);
v{d} = zeros(maxRank,1,p);
for j = 1:maxRank
    v{1}(1,j,:) = double(u{j}{1});
    v{d}(j,1,:) = double(u{j}{d});
end
for i = 2:d-1
    v{i} = zeros(maxRank,maxRank,p);
    for j = 1:maxRank
        v{i}(j,j,:) = double(u{j}{i});
    end
end


function [lopt,delta]=get_optimal_l(PC,fstemp,rstemp,tol,maxiter)
delta=zeros(length(fstemp),1);
dim=getnbgroups(PC);
Hs=polyval(PC,rstemp);
l=normalize(normalizephi(ones(1,1,PC)));
ls=cell(1,dim);
for i=1:dim
    ls{i} = Hs{i}*l{i};
end
for pp=1:30
    lp=l;
    for i=1:dim       
        l{0}=1;  
        ls{i}=ones(length(fstemp),1);
        gammas = prod([ls{:}],2);
        D = bsxfun(@times, Hs{i}, gammas);
        param.lambda=1.0e-03;
        param.numThreads=-1; 
        path=[];     
        [sol path]=mexLasso(full(fstemp),full(D),param);
        if nnz(path)>0
            [sol,dummy1,delta]=sel_opt_path_stat_test(D,fstemp,path);
            %[sol,dummy1,delta]=select_opt_path(D,fstemp,path,'leaveout','correction');
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
    if(error<1.0e-03)
        break;
    end
end
lopt=l;
