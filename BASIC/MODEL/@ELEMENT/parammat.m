function [se,fieldtype,fieldstorage,fieldddl] = parammat(elem,node,paramname,varargin)
% function [se,fieldtype,fieldstorage,fieldddl] = parammat(elem,node,paramname,varargin)

fieldtype = 'scalar';
fieldddl = DDL(paramname);
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
mat = getmaterial(elem);
se = evalparam(elem,mat,paramname,xnode,gauss.coord);

if ischarin('smooth',varargin)
    se = smooth(elem,node,se,gauss,'ddl',fieldddl);
end
