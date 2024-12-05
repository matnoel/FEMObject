function [u,v]=samesize(u,v)

n1= size(u,1);
n2= size(v,1);
if n1==1
u = repmat(u,n2,1);
end
if n2==1
v = repmat(v,n1,1);
end
if n1>1 & n2>1 & n1~=n2
   error('les points doivent etre de la meme taille') 
end

