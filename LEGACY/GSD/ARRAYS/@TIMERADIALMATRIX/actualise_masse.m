function A = actualise_masse(A,varargin)

if isempty(A.DLmasse) | length(A.DLmasse)~=A.m
A = calc_masse(A,varargin{:});
end
