function [se,fieldtype,fieldstorage,fieldddl] = opmatnode(elem,node,varargin)
% function [se,fieldtype,fieldstorage,fieldddl] = opmatnode(elem,node,varargin) 

[se,fieldtype,fieldstorage,fieldddl] = opmat(elem,node,'node',varargin{:});

se = transfernode(elem,node,se);
