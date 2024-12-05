function d=norm_residu(A,u,b,opt)
% function fastresdist(A,u,b)
%   || A*u-b ||/||b||
% Si opt=0 :
%   || A*u-b ||
if nargin==3
    opt=1;
end
if opt==1
    nb=norm(b);
else
    nb=1;
end

% stocker les alpha
ua=u.alpha;
Aa=A.alpha;
ba=b.alpha;
aualpha = ua'*Aa;
aualpha = [aualpha(:)' -ba];
% alpha=1
u.alpha=1;
u=gathervectors(u);
A.alpha=ones(size(Aa));
b.alpha=ones(size(ba));
% matrice de metrique
au=gathervectors(A*u + b);
au.F=cellfun(@(u) u'*u,au.F,'UniformOutput',0);
au=timesblock(au);
% erreur
d=sqrt(abs( aualpha*au*aualpha' ))/ nb;


A.alpha=Aa;
b.alpha=ba;