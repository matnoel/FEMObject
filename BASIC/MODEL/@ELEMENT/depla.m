function [se,fieldtype,fieldstorage,fieldddl] = depla(elem,node,q,varargin)
% function [se,fieldtype,fieldstorage,fieldddl] = depla(elem,node,q,varargin)

fieldtype = 'ddlnode';
fieldddl = getddlnode(elem);
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
N = calc_N(elem,xnode,gauss.coord);
se = N*qe;

if strcmp(getlstype(elem),'out')
    se(:)=0;
end
if ischarin('smooth',varargin)
    se = smooth(elem,node,se,gauss,'ddl',fieldddl);
end
