function hc = dpolyval(H,indices,j,x)
% function hc = dpolyval(H,indices,j,x)
% M : nombre de polynomes
% x : double n-by-M (n : nombre de valeurs)
% ou x : 1-by-M cell de n-by-1 double

M=getM(H);
indices = indices(:,1:M);
if isa(x,'cell')
    if length(x)~=M
        error('x doit etre une cell de length M')
    end
    x=[x{:}];
end
if size(x,2)~=M
    error('x doit etre avoir M colonne')
end

for k=setdiff(1:M,j)
    hc{k}=polyval(H.h{k},indices(:,k),x(:,k));
end
for k=j
    hc{k}=dpolyval(H.h{k},indices(:,k),x(:,k));
end

