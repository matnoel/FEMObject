function u = solve_updatenonleafnode_naive(A,b,u,t)
% function [u] = solve_updatenonleafnode_naive(A,b,u,t)
% Naive version, just used for benchmarking

N=numel(u.B{t});
Bz=zeros(size(u.B{t}));

uAu=zeros(N);
ub=zeros(N,1);
tic
for n=1:N
    fprintf('%2.2d \n',n/N);
    Br=Bz;
    Br(n)=1;
    ur=u;
    ur.B{t}=Br;
    Au=orthog(apply_mat_to_vec(A,ur));
    for m=1:N
        Bl=Bz;
        Bl(m)=1;
        ul=u;
        ul.B{t}=Bl;
        uAu(n,m)=innerprod(Au,ul);
    end
    ub(n)=innerprod(ur,b);
end
B=uAu\ub;

u.B{t}=reshape(B,size(u.B{t}));
end

