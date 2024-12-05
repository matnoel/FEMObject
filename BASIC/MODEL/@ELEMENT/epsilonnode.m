function [se,fieldtype,fieldstorage,fieldddl] = epsilonnode(elem,node,q,varargin)
% function [se,fieldtype,fieldstorage,fieldddl] = epsilonnode(elem,node,q,varargin) 

[se,fieldtype,fieldstorage,fieldddl] = epsilon(elem,node,q,'node',varargin{:});

se = transfernode(elem,node,se);
