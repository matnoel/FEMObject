function [apc,H]=decompfe(a,n,p,varargin)
% function [apc,H]=decompfe(a,n,p,'tol',tol)
% decomposition de a sur une base element fini
% n : nombre de subdivision 
% p : degre des polynomes
% tol : erreur souhaitee


if length(p)==1
   p=repmat(p,1,a.M);
end
if isa(n,'RANDPOLYS')
  n =  getpoly(n);
end
if ~isa(n,'cell')
n={n};
end
if length(n)==1 
   n=repmat(n,1,a.M);
end

for k=1:a.M
    [apc{k},h{k}] = decompfe(a.RV{k},'femesh',n{k},'order',p(k),varargin{:});       
end

H=RANDPOLYS(h);        
