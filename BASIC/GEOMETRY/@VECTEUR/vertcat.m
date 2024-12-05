function w = vertcat(u,v,varargin)
% function w = vertcat(u,v,varargin)

if isa(u,'VECTEUR') && isa(v,'VECTEUR')
    w = u;
    w.MYDOUBLEND = concat(u.MYDOUBLEND,v.MYDOUBLEND,3);
end

if length(varargin)>0
    w = vertcat(w,varargin{:});
end
