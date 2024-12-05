function varargout = gmsh2femobject(indim,file,varargin)
% function varargout = gmsh2femobject(indim,file,dim,option)
% indim : dimension de l'espace de travail
% file : nom du fichier .msh
% dim : recuperation des elements de dimension dim
%             (peut etre un vecteur)
%
% option :
% 'subfaces' :  si on veut que les facets, ridges, vertices aient leurs
% propres facets, ridges, vertices.
% 'nofaces' : si on ne veut pas de facets, ridges, peaks (pour sauver de l'espace memoire)

file = GMSHFILE(file);
n = max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,file,varargin{:});
