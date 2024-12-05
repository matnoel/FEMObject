function [taue,fieldtype,fieldstorage,fieldddl] = stabfactor(elem,node,varargin)

fieldtype = 'scalar';
fieldstorage = 'center';
fieldddl = DDL('tau');

gauss=calc_gauss(elem,0);
xnode = node(elem);
mat=getmaterial(elem);

taue = stabfactor(mat,elem,xnode,gauss.coord);

