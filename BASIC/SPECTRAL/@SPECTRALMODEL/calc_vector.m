function b = calc_vector(S,varargin)
% function b = calc_vector(S,varargin)

nofree = ischarin('nofreematrix',varargin);
varargin = delonlycharin('nofreematrix',varargin);

if S.dim==1
    b = calc_vector(S.L,varargin{:});
elseif S.dim==2
    b = calc_vector2D(S,varargin{:});
end
if ~nofree
    b = freevector(S,b);
end
