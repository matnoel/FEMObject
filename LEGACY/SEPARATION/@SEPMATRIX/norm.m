function n = norm(u,varargin)
% function n = norm(u,varargin)

if isa(u,'SEPMATRIX')
    if nargin == 1
        alpha = u.alpha;
        u.alpha = 1;
        u = gathervectors(u);
        u = timesblock(u'*u);
        n = sqrt(abs(alpha*u*alpha'));
    elseif nargin == 2
        alpha = u.alpha;
        u.alpha = 1;
        u = gathervectors(u);
        u = timesblock(u'*varargin{1}*u);
        n = sqrt(abs(alpha*u*alpha'));
    end
elseif isa(u,'double') && nargin==2 && isa(varargin{1},'double')
    n = sqrt(u'*varargin{1}*u);
end

end


function n = norm_old(u,varargin)

if nargin==2
    
    if isa(u,'SEPMATRIX')
        v = varargin{1}*u;
    elseif isa(u,'double')
        Ametric = varargin{1};
        if Ametric.m>1
            error('ne marche qu''avec un rang de 1 pour la matrice metrique')
        end
        v = Ametric.alpha*u;
        dim = getdim(Ametric);
        for j=1:length(Ametric.F)
            rep=1:dim;
            rep(1)=j;
            rep(j)=1;
            v = permute(v,rep);
            v = Ametric.F{j}*v;
            v = permute(v,rep);
        end
    else
        error('pas prevu');
    end
else
    v = u;
end


if isa(u,'double')
    n = full(abs(sqrt(prodscal(u,v))));
else
    if size(u.F{1},2)>1
        v = times(u,v);
        v = sum(sum(v));
    else
        v = mtimes(u',v);
    end
    
    n = full(reshape([v.F{:}],v.m,v.dim));
    n = sqrt(abs(sum(prod(n,2).*v.alpha(:))));
end

end