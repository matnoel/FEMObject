function [optpattern, output]=sel_opt_pattern_stat_test(D,y,pattern)
% function [optpattern, output]=sel_opt_pattern_stat_test(D,y,pattern)
% Selection of a sparsity pattern using cross-validation (loo)
% D is N-by-P, y is of size N-by-1
% pattern : P-by-m matrix whose columns give m different sparsity patterns 
% (is pattern is not boolean, then pattern <- (pattern~=0))
% optpattern: a boolean vector indicating the optimal sparsity pattern 
% output.relativeError : corrected loo relative error estimator
% output.ind : 

N=numel(y);

pattern(:,1)=[];
X=(pattern~=0);
[~, I]=unique(X','rows','first');
pattern=pattern(:,sort(I));


xpattern=zeros(size(pattern));
for i=1:size(pattern,2)
    sol=pattern(:,i);
    ind=find(sol);
    Dred=D(:,ind);
    x=(Dred'*Dred)\(Dred'*y);
    vec=zeros(size(sol,1),1);
    vec(ind)=x;
    xpattern(:,i)=vec;
end


epsil=zeros(size(pattern,2),1);

GG = D'*D;
for i = 1:size(xpattern,2)
    Dx=D*xpattern(:,i);
    ind=find(pattern(:,i));
    P=D(:,ind);
    G=inv(GG(ind,ind));    
    C = sum(P'.*(G*P'),1)';
    merror=((y-Dx)./(1-C)).^2;
    error=sum(merror,1)/N;
    epsil(i) = error/var(y);
    corrFact=(1/(1-nnz(xpattern(:,i))/N))*(1+trace(G));
    epsil(i)=epsil(i)*corrFact;
end
[~, ind] = min(epsil);
optpattern = pattern(:,ind);

output.relativeErrorpattern = sqrt(epsil);
output.relativeError = sqrt(epsil(ind));
output.ind = ind-1;







