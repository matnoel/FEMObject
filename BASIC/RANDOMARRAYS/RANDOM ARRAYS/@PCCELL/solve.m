function [u,varargout] = solve(A,b,varargin)
% function [u,varargout] = solve(A,b)
% Resolution de Ax=b
% A : PCCELL
% b : PCARRAY ou PCCELL ou PCRADIAL
%
% function [u,varargout] = solve(A,b,PC)
% PC : POLYCHAOS sur lequel est defini b et la solution x (A peut etre defini sur un autre)
%
% function [u,varargout] = solve(A,b,PC,'norandom',norandom)
% les composantes norandom de x sont deterministes


nosto = getcharin('norandom',varargin);
PC = getclassin('POLYCHAOS',varargin,getPC(A));
PC = calc_masse(PC,getPC(A));
P = getP(PC);

if isa(b,'double')
    b = b.*one(PC);
end
b = double(expand(b));
n = size(b,1);
nns = length(nosto);
sto = setdiff(1:n,nosto);
ns = length(sto);
b1 = b(sto,:);
b2 = b(nosto,:)*double(one(PC))';

masse = getmasse(PC);
for k=1:length(masse)
    masseun{k} = masse{k}*double(one(PC))';
end
massep = getmatrix(masse,1:P+1,1:P+1);
masseunp = getmatrix(masseun,1:P+1,1);
val11 = getmatrix(A.value,sto,sto);
val12 = getmatrix(A.value,sto,nosto);
val21 = getmatrix(A.value,nosto,sto);
val22 = getmatrix(A.value,nosto,nosto);


val11 = val11*massep';
val12 = val12*masseunp';
val21 = val21*masseunp';
val22 = val22*double(one(getPC(A)))';

M11 = assemble(val11,ns,ns,P+1,P+1);
M12 = assemble(val12,ns,nns,P+1,1);
M21 = assemble(val21,nns,ns,1,P+1);
M22 = assemble(val22,nns,nns,1,1);

u = zeros(n,P+1);
v = [M11,M12;M21,M22]\[b1(:);b2(:)];
v1 = v(1:ns*(P+1));
v2 = v(ns*(P+1)+1:end);


u(sto,:) = reshape(v1,ns,P+1);
u(nosto,:) = v2*double(one(PC));

u = PCARRAY(u,PC);



function M = assemble(val,n1,n2,p1,p2)
% function M = assemble(val,n1,n2,p1,p2)

[I,J] = find(val);
[Ji,Jj] = ind2sub([p1,p2],J);
[Ii,Ij] = ind2sub([n1,n2],I);
j = (Jj(:)-1)*n2+Ij(:);
i = (Ji(:)-1)*n1+Ii(:);
I = sub2ind(size(val),I,J);
M = sparse(i,j,val(I),n1*p1,n2*p2);

return


