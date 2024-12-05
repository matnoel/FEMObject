function [se,fieldtype,fieldstorage,fieldddl] = sigmanode(elem,node,q,varargin)
% function [se,fieldtype,fieldstorage,fieldddl] = sigmanode(elem,node,q,varargin)

[se,fieldtype,fieldstorage,fieldddl] = sigma(elem,node,q,'node',varargin{:});

se = transfernode(elem,node,se);
