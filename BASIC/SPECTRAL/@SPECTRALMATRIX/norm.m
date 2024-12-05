function an = norm(a,varargin)
% function an = norm(a,varargin)

an = full(sqrt(abs(prodscal(a,a,varargin{:}))));

% if nargin==2
%     if ~israndom(A)
%         b = A*a.MULTIMATRIX;
%     else
%         b = A*a;
%         
%         b = b.MULTIMATRIX;
%     end
% else
%     b=a.MULTIMATRIX;
% end

% an = b(:)'*a.MULTIMATRIX(:);
% 
% an = full(sqrt(sum(double(an))));
