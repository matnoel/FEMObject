function x = MYDOUBLEND(u,k)
% function x = MYDOUBLEND(u)
% u : MULTIMATRIX
% x : MYDOUBLEND
%
% function x = MYDOUBLEND(u,k)
% k : indique la dimension dans x correspondant aux multimatrix

if nargin==1
 x = MYDOUBLEND(reshape(full(u.value),[u.s,u.sm])); 
else
 s = u.s ;
 if length(k)==1
 sm = prod(u.sm);
 else
 sm = u.sm ;   
 end
 x = MYDOUBLEND(reshape(full(u.value),[s,sm])); 
 
 p=[1:max(length(s)+length(sm),max(k))];
 locm = length(s)+[1:length(sm)];
 p(k)=locm;
 p(locm)=k;
 
 x = permute(x,p);
 
 
end