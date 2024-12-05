function u = solve_alphaupdate_sep(A,b,W,Wtilde)
%function u = solve_alphaupdate_sep(A,b,W,Wtilde)

W.alpha=ones(1,W.m);
if nargin==3
    Wtilde = W;
else
    Wtilde.alpha=ones(1,Wtilde.m);
end

AW=cell(W.m,1);
Wb=zeros(W.m,1);

parfor i=1:W.m
    Wi=TSEPMATRIX(truncate(W,i));
    AW{i}=A*Wi;
    Wb(i)=expand(Wi'*b);
end

WAW=zeros(Wtilde.m,W.m);

for j=1:Wtilde.m
    parfor i=1:W.m
        WAW(j,i)=expand(...
            TSEPMATRIX(truncate(Wtilde,j))'*AW{i});
    end
end

lambda=WAW\Wb;

u=W;
u.alpha=lambda';
