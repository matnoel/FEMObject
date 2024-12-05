function [U,xie,True_error_ex]=greedy_exact(A,b,rank,Xi)


xie=zeros(rank+1,1);
True_error_ex=zeros(rank,1);

[~,xie(1)]=min(sum((Xi-0.5).^2,2));

for r=1:rank
    fprintf('%d ',r)
    % Nouveau vecteur dans la base:
    v= eval_sparse(A,xie(r))\b;
    if r==1
        U = v;
    else
        U = [U v];
    end
    % on orthogonalise pour le fun
    [U,~,~]=svd(U,'econ');
    % recherche du nouveau point
    parfor k=1:A.n
        uk = eval_sparse(A,k)\b;
        u_proj = U*(U'*uk);
        error(k)=norm(uk-u_proj)/norm(uk);
    end
    
    [True_error_ex(r),xie(r+1)]=max(error);
    
    
    
end

