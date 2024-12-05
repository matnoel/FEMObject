function w = horzcat(u,v,varargin)
% function w = horzcat(u,v,varargin)

if  isa(u,'VECTEUR') && isa(v,'VECTEUR')
    w = u;
    w.MYDOUBLEND = concat(u.MYDOUBLEND,v.MYDOUBLEND,4);
end

if length(varargin)>0
    w = horzcat(w,varargin{:});
end
