function [se,fieldtype,fieldstorage,fieldddl] = energyintnode(elem,node,q,varargin)
% function [se,fieldtype,fieldstorage,fieldddl] = energyintnode(elem,node,q,varargin) 

[se,fieldtype,fieldstorage,fieldddl] = energyint(elem,node,q,'node',varargin{:});

se = transfernode(elem,node,se);


