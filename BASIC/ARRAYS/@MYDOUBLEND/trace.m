function t = trace(A)
%TRACE  Sum of diagonal elements.
%   TRACE(A) is the sum of the diagonal elements of A, which is
%   also the sum of the eigenvalues of A.
%
%   Class support for input A:
%      MYDOUBLEND

if size(A,1)~=size(A,2)
    error('Matrix must be square.');
end
t = full(sum(diag(A)));
