function W = mygram(W,tol)
if nargin==1
    tol=1e-15;
end
q=0;
for i=1:size(W,2)
    w = W(:,i);
    w = w/sqrt(w'*w);
    w = w - W(:,1:q)*(W(:,1:q)'*w);
    n = sqrt(w'*w);
    if n>tol
        q=q+1;    
        w = w/n;
        W(:,q) = w;
%    warning('critical orthogonalization')
    end
end

W=W(:,1:q);
return
