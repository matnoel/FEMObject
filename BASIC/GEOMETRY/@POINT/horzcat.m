function w = horzcat(u,v,varargin)
% function w = horzcat(u,v,varargin)

if isa(u,'POINT') && isa(v,'POINT')
    w = u;
    w.MYDOUBLEND = concat(u.MYDOUBLEND,v.MYDOUBLEND,4);
end

if ~isempty(varargin)
    w = horzcat(w,varargin{:});
end