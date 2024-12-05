function u = remesh(u,options,varargin)
% function u = remesh(u,dim)
% launch gmsh u -dim

% function u = remesh(u,[],'gmshoptions',options)
% options: char
% launch gmsh u options

if nargin==1
    error('rentrer les options de gsmh')
end

if ~isempty(options) && isa(options,'double')
    options = [' -' num2str(options)];
else
    options =  '';
end

if ischarin('gmshoptions',varargin)
    options = [options ' ' getcharin('gmshoptions',varargin)];
end

if isempty(strfind(options,'-format')) % for compatibility with Matlab version < 9.1 (R2016b)
% if ~contains(options,'-format') % for Matlab versions >= 9.1 (R2016b)
    options = [options ' -format msh2'];
end

u = runfilemsh(u,options);
