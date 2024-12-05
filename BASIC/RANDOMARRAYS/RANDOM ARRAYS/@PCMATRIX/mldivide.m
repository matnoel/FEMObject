function [u,varargout] = mldivide(A,b)
varargout=cell(1,nargout-1);
[u,varargout{:}] = solve(A,b);