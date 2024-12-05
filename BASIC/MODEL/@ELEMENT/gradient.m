function [se,fieldtype,fieldstorage,fieldddl] = gradient(elem,node,q,varargin)
% function [se,fieldtype,fieldstorage,fieldddl] = gradient(elem,node,q,varargin)

fieldtype = 'ddlgauss';
fieldddl = getddlgauss(elem);
if ischarin('node',varargin)
    gauss.coord = permute(nodelocalcoord(elem),[4,2,3,1]);
    gauss.nbgauss = getnbnode(elem);
    fieldstorage = 'node';
elseif ischarin('smooth',varargin)
    n = getcharin('intorder',varargin,'mass');
    gauss = calc_gauss(elem,n);
    fieldstorage = 'node';
else
    n = getcharin('intorder',varargin,'rigi');
    gauss = calc_gauss(elem,n);
    fieldstorage = 'gauss';
end

xnode = getcoord(node,getconnec(elem)');

qe = localize(elem,q);
B = calc_gradient(elem,xnode,gauss.coord);
se = B*qe;

if ischarin('smooth',varargin)
    se = smooth(elem,node,se,gauss,'type',fieldtype);
end

