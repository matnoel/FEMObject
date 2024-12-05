function varargout = calc_mode(K,M,varargin)
% function [V,D] = calc_mode(K,M,varargin)

rep= K.ddlfree;

A = K.MYDOUBLE(rep,rep);
B = M.MYDOUBLE(rep,rep); B = (B+B')/2;
m=varargin{1};
varargout{:} = eigs(A,B,max(m),'SM');

if nargout ==1
    D = eigs(A,B,max(m),'SM');
    D = sort(D);
    varargout{1} = D(m);
elseif nargout==2
    [Vrep,D] = eigs(A,B,max(m),'SM');
    [D,renum] = sort(diag(D));
    D = diag(D(m));
    V = zeros(size(K,1),length(m));
    V(rep,:) = Vrep(:,renum(m));
    varargout{1} = FEVECTOR(V,K.ddlbloque);
    varargout{2} = D;
end
