function [se,fieldtype,fieldstorage,fieldddl] = parammatnode(elem,node,paramname,varargin)
% function [se,fieldtype,fieldstorage,fieldddl] = parammatnode(elem,node,paramname,varargin) 

[se,fieldtype,fieldstorage,fieldddl] = parammat(elem,node,paramname,'node',varargin{:});

se = transfernode(elem,node,se);
