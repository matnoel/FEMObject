function A = actualise_ximasse(A,varargin)

if isempty(A.DLmasse) || length(A.DLmasse)~=A.m 
A=calc_ximasse(A,varargin{:});
end
