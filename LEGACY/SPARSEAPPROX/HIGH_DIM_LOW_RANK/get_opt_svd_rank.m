function [rankopt]=get_opt_svd_rank(A,fs,s,maxrank,reg)

Kfold=3;
[fsReg,AReg,fsTest,ATest]=cross_validation_Kfold_least_squares(A,fs,Kfold);
error=zeros(maxrank,Kfold);
for k=1:Kfold
    param.lambda=1.0e-02;
    if (strcmp(reg,'ols')==1)
        solopt=(AReg{k}'*AReg{k})\(AReg{k}'*fsReg{k});
    else
        [~,path]=mexLasso(full(fsReg{k}),full(AReg{k}),param);
        solopt=sel_opt_path_stat_test(AReg{k},fsReg{k},path);
        %[solopt,epsil]=select_opt_path(AReg{k},fsReg{k},path,'leaveout');
    end
    W=reshape(solopt,s(1)*s(2),s(3)*s(4));
    [U,Si,V]=svd(full(W));
    
    
    for tt=1:maxrank
        if tt>numel(diag(Si))
            break
        end
        S=U(:,1:tt)*Si(1:tt,1:tt)*V(:,1:tt)';
        fapproxpath=ATest{k}*reshape(S,size(S,1)*size(S,2),1);
        error(tt,k)=norm((fsTest{k}-fapproxpath),2)/norm(fsTest{k},2);
    end
end
    [~,rankopt]=min(sum(error,2)/Kfold);
   
    