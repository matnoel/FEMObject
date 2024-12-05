function hc = polyval(H,indices,x)
% function hc = polyval(H,indices,x)
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
    error('x doit etre avoir M colonnes')
end

for k=1:M
    hc{k}=polyval(H.h{k},indices(:,k),x(:,k));
end

