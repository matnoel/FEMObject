function [u,varargout] = solve_singular(A,b,varargin)
% function [u,varargout] = solve(A,b)
% Resolution de Ax=b
% A : PCMATRIX
% b : PCMATRIX ou PCRADIALMATRIX ou double
% 
% function [u,varargout] = solve(A,b,PC)
% PC : POLYCHAOS sur lequel est defini b et la solution x (A peut etre defini sur un autre)
%
% function [u,varargout] = solve(A,b,PC,'norandom',norandom)
% les composantes norandom de x sont deterministes



[A,b,PC]=pcsystemupdate(A,b,varargin{:});


n=size(b,1);
Amean = expect(A);
repzero = find(abs(diag(Amean))/max(abs(diag(Amean)))<1e-12);
repnonzero = setdiff(1:n,repzero);
varargout=cell(1,nargout);

if ~isempty(repzero)
    [utemp,varargout{:}]=solve(getcompo(A,repnonzero,repnonzero),getcompo(b,repnonzero,1),varargin{:});
     u = zeros(n,length(PC));
 u(repnonzero,:)=double(utemp);
 u = PCMATRIX(u,[n,1],PC);
else
    [u,varargout{:}] = solve(A,b,varargin{:});
end

    