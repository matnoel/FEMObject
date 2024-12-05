function [u] = solve_updatenonleafnode(A,b,u,t)
% function [u] = solve_updatenonleafnode(A,b,u,t)

N=numel(u.B{t});
uAu=zeros(N);
Au=apply_mat_to_vec(A,u);

sz_A = size(A.B{t});
sz_A(end+1:3) = 1;
sz_u = size(u.B{t});
sz_u(end+1:3) = 1;

M = reduceinnerprodatnode(u,Au,t);
Mb = reduceinnerprodatnode(u,b,t);

if t==1 %t is root node
    for n=1:N
        Br=zeros(size(u.B{1}));
        Br(n)=1;
        AuB=Au.B{1};
        AuB(:, :) = kron(A.B{1}(:, :), Br(:, :));
        B = M{1}*AuB*M{2}';
        uAu(:,n)=B(:);
    end
    ub = ttm(b.B{1}, Mb(1:2), 1:2);
    ub=ub(:);
else
    Bb=ttm(b.B{t},Mb,[1 2 3]);
    szuBt=size(u.B{t}); %Suppressing "." make things faster
    szAuB=size(Au.B{t});
    AuB=zeros(szAuB);
    iA=((1:sz_A(1)) - 1) * sz_u(1);
    jA=((1:sz_A(2)) - 1) * sz_u(2);
    kA=((1:sz_A(3)) - 1) * sz_u(3);
    ABt=A.B{t};
    for n=1:N
        if n>1
% Allocating at each iteration takes too much time, 
% almost as much as "B = ttm(AuB,M,[1 2 3])".
% Erasing previous entries is much faster.
            AuB(iAuB,jAuB,kAuB)=zeros(numel(iAuB),numel(jAuB),numel(kAuB));
        end
        [iu,ju,ku] = ind2sub(szuBt,n);
        iAuB=iu+iA;
        jAuB=ju+jA;
        kAuB=ku+kA;
        AuB(iAuB,jAuB,kAuB)=ABt;
        B = ttm(AuB,M,[1 2 3]);
        uAu(n,:) = B(:);
    end
    ub = Bb(:);
end

u.B{t} = reshape(uAu\ub,size(u.B{t}));

end

