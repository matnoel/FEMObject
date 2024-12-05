function [T,err] = hooi(X,r,varargin)
% function T = hooi(X,r,varargin)
% Higher-order orthogonal iteration
% X : N-way tensor
% r : rank of the HOSVD
% T : TSEPMATRIX
% Tensor decompositions and applications - Kolda, Bader, SIAM, 2009

N=ndims(X);
sz=size(X);

if length(r)==1
    r=r*ones(1,N);
elseif length(r)~=N
    error('not implemented')
end

if ischarin('SEPSOLVER',varargin)
    solver = getcharin('SEPSOLVER',varargin);
else
    solver = SEPSOLVER(N,varargin{:});
end
param = getparam(solver);

if param.maxiter==10
    param.maxiter=50;
end

erriter = zeros(1,param.maxiter+1);
[t,err]=hosvd(X,r);
erriter(1)=err;
if err>param.tol
    gt=gathervectors(t);
    A=gt.F;
    
    normX=norm(X(:));
    
    for k=1:param.maxiter
        for n=1:N
            dims=[1:n-1,n+1:N];
            Yn=ttm(X,A{dims(1)}',dims(1));
            for i=2:N-1
                Yn=ttm(Yn,A{dims(i)}',dims(i));
            end
            Yn=permute(Yn,[n dims]);
            Yn=reshape(double(Yn),sz(n),prod(r(dims)));
            [A{n},temp1,temp2]=svds(Yn,r(n));
        end
        alpha=ttm(X,A{1}',1);
        for n=2:N
            alpha=ttm(alpha,A{n}',n);
        end
        erriter(k+1)=(normX^2-norm(alpha(:))^2)/normX^2;
        stag=abs(erriter(k+1)-erriter(k))/erriter(k+1);
        if param.display
            fprintf('  iteration #%d - residual = %d - stagnation = %d\n',k,erriter(k+1),stag)
        end
        if stag<param.tol
            break
        end
    end
    
    err=abs((normX^2-norm(alpha(:))^2))/normX^2;
    if param.display
        fprintf('  final error = %d\n',err)
    end
    
    T=TSEPMATRIX(N);
    T.alpha=alpha;
    
    for n=1:N
        for j=1:r(n)
            T.F{n}{j}=A{n}(:,j);
        end
    end
else
    T=t;
end
