function w = norm(u,normtype)
% function w = norm(u,normtype)
% norme du vecteur u
% normtype : 1, 2 ou Inf (type de norme) 2 par defaut
%
% See also VECTEUR/dot, VECTEUR/norm, VECTEUR/normalize, VECTEUR/cross, VECTEUR/abs, VECTEUR/planortho, VECTEUR/minus, VECTEUR/rot2D,
% VECTEUR/mtimes, VECTEUR/times, VECTEUR/mrdivide, VECTEUR/ne, VECTEUR/eq,
% VECTEUR/plus, VECTEUR/uminus, VECTEUR/norm, POINT/norm

if nargin==1
    normtype = 2;
end

switch normtype
    case 2
        w = sqrt(sum(u.MYDOUBLEND.^2,1));
    case 1
        w = sum(abs(u.MYDOUBLEND),1);
    case Inf
        w = max(abs(u.MYDOUBLEND),[],1);
end
