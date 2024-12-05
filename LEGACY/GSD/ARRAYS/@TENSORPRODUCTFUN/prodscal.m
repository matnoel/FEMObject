function s = prodscal(a,b,V)
% produit scalaire de deux fonctions TENSORPRODUCTFUN a et b
% d : produit scalaire suivant la dimension d


if nargin==2 && a.PRODUCTSPACE~=b.PRODUCTSPACE
error('indiquer les dimensions pour le produit scalaire, sinon, les TENSORPRODUCTFUN doivent vivre sur le meme espace produit')
elseif nargin==2 && a.PRODUCTSPACE==b.PRODUCTSPACE
 V = a.PRODUCTSPACE;
elseif nargin==3 && ~isa(V,'PRODUCTSPACE')
 V = PRODUCTSPACE(V);
end
Vglob = union(a.PRODUCTSPACE,b.PRODUCTSPACE);     

if ~all(ismember(V,Vglob))
    error('les dimensions indiquees n''existent pas dans les TENSORPRODUCTFUN')
end
[sglob,repglob] = ismember(Vglob,V);

s.factor = a.factor*b.factor;
s.phi = cell(1,getnbdim(Vglob));
for i=1:getnbdim(Vglob)
[sa,repa] = ismember(Vglob(i),a.PRODUCTSPACE);
[sb,repb] = ismember(Vglob(i),b.PRODUCTSPACE);

if sglob(i)               
if sa && sb 
s.factor = s.factor * prodscal(a.phi{repa},b.phi{repb});
elseif sa
s.factor = s.factor * prodscal(a.phi{repa},1); 
elseif sb
s.factor = s.factor * prodscal(1,b.phi{repb}); 
else
    error('pas prevu')
end
else
if sa && sb 
s.phi{i} = mtimes(a.phi{repa},b.phi{repb});
elseif sa
s.phi{i} = a.phi{repa}; 
elseif sb
s.phi{i} = b.phi{repb}; 
else
    error('pas prevu')
end    
    
end
end
s.phi(sglob)=[];

s = TENSORPRODUCTFUN(setdiff(Vglob,V),s.factor,s.phi{:});
s = simplify(s);



