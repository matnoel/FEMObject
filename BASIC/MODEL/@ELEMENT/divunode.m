function [se,fieldtype,fieldstorage,fieldddl] = divunode(elem,node,q,varargin)
% function [se,fieldtype,fieldstorage] = divunode(elem,node,q,varargin) 

[se,fieldtype,fieldstorage,fieldddl] = divu(elem,node,q,'node',varargin{:});

se = transfernode(elem,node,se);
