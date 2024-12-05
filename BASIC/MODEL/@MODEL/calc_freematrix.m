function K=calc_freematrix(S,varargin)

if israndom(S)
    K = calc_pcmatrix(S,varargin{:});
else
    K = calc_matrix(S,varargin{:});
end
K = freematrix(S,K);