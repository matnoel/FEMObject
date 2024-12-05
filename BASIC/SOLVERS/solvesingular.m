function [x,flag]=solvesingular(A,b,tol,varargin)
% function [x,flag]=solvesingular(A,b,tol,varargin)

if nargin<=2 || isempty(tol)
    tol=1e-12;
end
d=diag(A);
rep = find(abs(d)./max(abs(d))>tol);
x=zeros(size(b));

x(rep,:)=A(rep,rep)\b(rep,:);
flag=0;
% fprintf('condest')
% condest(A(rep,rep))
