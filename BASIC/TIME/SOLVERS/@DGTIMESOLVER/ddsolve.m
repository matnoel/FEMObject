function [ut,result,vt,at] = ddsolve(L,b,M,K,C,u0,v0,varargin)
% function [ut,result,vt,at] = ddsolve(L,b,M,K,C,u0,v0,solver)
% Discontinuous Galerkin solver: solves Mu''+Cu'+Ku=b
% L: DGTIMESOLVER
% b: right-hand side
% M,K,C: double or random matrices
% solver: function for solving linear systems
%         solver(A,b) solves Au=b
% ut: TIMEMATRIX of solution u
% vt: TIMEMATRIX of velocity v=u'
% at: TIMEMATRIX of acceleration a=v'=u''
% result: struct containing outputs
% result.totaltime : total CPU time

% display_ = getparam(L,'display');

n = size(M,1);
if nargin<5 || isempty(C)
    C = sparse(n,n);
end

[b,M,K,C] = init_resolution(L,b,M,K,C);

A = [M,C;sparse(n,n),M];
B = [sparse(n,n),K;-M,sparse(n,n)];
b0 = [b;(sparse(n,1)*one(L))];

if (nargin<6 || isempty(u0)) && (nargin<7 || isempty(v0))
    u0 = [];
else
    if nargin<6 || isempty(u0)
        u0 = sparse(n,1);
    end
    if nargin<7 || isempty(v0)
        v0 = sparse(n,1);
    end
    u0 = [v0;u0];
end

[Ut,result] = dsolve(L,b0,A,B,u0,varargin{:});

if nargout>=4
    % Split ut and vt
    vt = Ut(1:n);
    ut = Ut(n+1:2*n);
    % Compute at
    at = dderivative(L,ut,vt,b,M,K,C,varargin{:});
end

if ~getparam(L,'outputsplit')
    % Do not split ut and vt
    ut = Ut;
elseif nargout<4
	% Split ut and vt if it isn't done yet
    vt = Ut(1:n);
    ut = Ut(n+1:2*n);
end


% % Change of variables
% % Compute ut, vt and at simultaneously with a change of variables
% % This method is slower than post processing.
% % To use this method, add the a0 input argument to the function definition:
% % function [ut,result,vt,at] = ddsolve(L,b,M,K,C,u0,v0,a0,varargin)
% 
% n = size(M,1);
% if nargin<5 || isempty(C)
%     C = sparse(n,n);
% end
% 
% [b,M,K,C] = init_resolution(L,b,M,K,C);
% 
% A = [C,M,sparse(n,n);M,sparse(n,n),sparse(n,n);sparse(n,n),M,sparse(n,n)];
% B = [K,sparse(n,n),sparse(n,n);sparse(n,n),-M,sparse(n,n);sparse(n,n),sparse(n,n),-M];
% b = [b;(sparse(2*n,1)*one(L))];
% 
% if (nargin<6 || isempty(u0)) && (nargin<7 || isempty(v0)) && (nargin<8 || isempty(a0))
%     u0 = [];
% else
%     if nargin<6 || isempty(u0)
%         u0 = sparse(n,1);
%     end
%     if nargin<7 || isempty(v0)
%         v0 = sparse(n,1);
%     end
%     if nargin<8 || isempty(a0)
%         a0 = sparse(n,1);
%     end
%     u0 = [v0;u0;a0];
% end
% 
% [Ut,result] = dsolve(L,b,A,B,u0,varargin{:});
% 
% if ~getparam(L,'outputsplit')
%     % Do not split ut and vt
%     ut = Ut;
% else
%     % Split ut, vt and at
%     ut = Ut(1:n);
%     vt = Ut(n+1:2*n);
%     at = Ut(2*n+1:3*n);
% end
