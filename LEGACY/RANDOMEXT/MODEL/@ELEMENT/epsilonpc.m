function [se,fieldtype,fieldstorage,fieldddl] = epsilonpc(elem,node,q,varargin)
fieldtype = 'ddlgauss';
fieldddl=getddlgauss(elem);
if ischarin('node',varargin)
gauss.coord=permute(nodelocalcoord(elem),[4,2,3,1]);
gauss.nbgauss = getnbnode(elem);
fieldstorage = 'node';
elseif ischarin('smooth',varargin)
n=getcharin('intorder',varargin,'mass');
gauss=calc_gauss(elem,n);        
fieldstorage = 'node';
else
    
n=getcharin('intorder',varargin,'rigi');
gauss=calc_gauss(elem,n);        
fieldstorage = 'gauss';

end

mat=getmaterial(elem);

xnode = getcoord(node,getconnec(elem)');

[qe,L]=localizepc(elem,q); 

[B,detJ]=calc_B(elem,xnode,gauss.coord);
se=B*qe;


if ischarin('smooth',varargin)
if size(se,2)~=1
    error('pas programme')
end
if israndom(q)
se = permute(se,[1,5,3,4,2]);
se = smooth(elem,node,se,gauss);
se = permute(se,[1,5,3,4,2]);
else
se = smooth(elem,node,se,gauss);
end
end
 
 if isa(se,'PCMYDOUBLEND')
    return 
 elseif  isa(L,'POLYCHAOS')
 se = PCMATRIX(MULTIMATRIX(se,5),L);
 elseif isa(L,'PCMATRIX')
 se = MULTIMATRIX(se,5);
 se = PCRADIALMATRIX(se,size(se),L);    
 end

