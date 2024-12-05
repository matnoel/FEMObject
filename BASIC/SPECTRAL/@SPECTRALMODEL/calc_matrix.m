function B = calc_matrix(S,varargin)
% function B = calc_matrix(S,varargin)

nofree = ischarin('nofreematrix',varargin);
varargin = delonlycharin('nofreematrix',varargin);

if S.dim==1
    B = calc_matrix(S.L,varargin{:});
else
    B = calc_matrix2D(S,varargin{:});
end

if isa(S.L,'POLYFELAGRANGE')
    B = sparse(B);
end

if ~nofree
    B = freematrix(S,B);
end


