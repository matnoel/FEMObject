function [U,maxNPr,xie,ResNorm,True_error]=greedy_preconditioned_RB(A,b,P,m,Xi)



xie=zeros(m+1,1);


[~,xie(1)]=min(sum((Xi-0.5).^2,2));



maxNPr =zeros(m,1);
ResNorm = cell(m,1);
True_error = cell(m,1);

for r=1:m
    fprintf('%d ',r)
    % Nouveau vecteur dans la base:
    v = eval_sparse(A,xie(r))\b;
    if r==1
        U = v;
    else
        U = [U v];
    end
    % on orthogonalise pour le fun
    [U,~,~]=svd(U,'econ');
    % preparatifs
    PAU = A*U;
    PAU = mtimes_truncation(P,PAU);
    PB  = P*b;
    residual=zeros(A.n,1);
    error=zeros(A.n,1);
    % recherche du nouveau point
    parfor k=1:A.n
        PAU_k = eval(PAU,k);
        PB_k  = eval(PB,k);
        u_r = (U'*PAU_k) \ (U'*PB_k);
        % norme du residu
        residual(k) = norm( PAU_k*u_r -PB_k)/norm(PB_k) ;
        
        % erreur vraie
%         if nargout==5
            uk = eval_sparse(A,k)\b;
            error(k)=norm(uk-U*u_r)/norm(uk);
%         end
        
    end
    ResNorm{r}=residual;
    True_error{r}= error;
    
    
    [maxNPr(r),xie(r+1)]=max(residual);
    
end
disp(' ')





