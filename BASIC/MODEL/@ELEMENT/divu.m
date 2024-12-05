function [se,fieldtype,fieldstorage,fieldddl] = divu(elem,node,q,varargin)
% function [se,fieldtype,fieldstorage,fieldddl] = divu(elem,node,q,varargin)

fieldtype = 'scalar';
fieldddl = DDL('DIVU');
if ischarin('node',varargin)
    gauss.coord = permute(nodelocalcoord(elem),[4,2,3,1]);
    gauss.nbgauss = getnbnode(elem);
    fieldstorage = 'node';
elseif ischarin('smooth',varargin)
    n = getcharin('intorder',varargin,'mass');
    gauss = calc_gauss(elem,n);
    fieldstorage = 'node';
else
    gauss = calc_gauss(elem,1);
    fieldstorage = 'gauss';
end

xnode = getcoord(node,getconnec(elem)');

qe = localize(elem,q);
DN = calc_DN(elem,xnode,gauss.coord);
se = DN(1,:)*qe(1:2:end,:)+DN(2,:)*qe(2:2:end,:);

if ischarin('smooth',varargin)
    se = smooth(elem,node,se,gauss,'ddl',fieldddl);
end
