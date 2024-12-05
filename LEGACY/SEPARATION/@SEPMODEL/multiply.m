function w=multiply(SM,u,v,SMw,varargin)
% Realise la des fonctions w=u*v a partir de la connaissance de SM
% ATTENTION : u et v du meme format


if isa(u,'HSEPMATRIX')
    % w=masse*u*v
    w=mtimeslike(HSEP(SM.masse ,u.tree),u,@SEPtimesSEP_FH_MULTM,@HSEPMATRIX)*v;
    % w=metric\w
    if SM.metric.m==1
        w=mtimeslike(HSEP(SM.metric,w.tree),w,@mldivide);
    elseif nargin>3 && isa(varargin,'HSEPSOLVER')
        w=solve(HSEP(SM.metric,w.tree),varargin{1});
    else
        error('Donner un HSEPSOLVER en argument')
    end
    % Demande de recompression de donnees
    if nargin>3, w=multisvd(w,varargin{:});end
elseif isa(u,'SEPMATRIX')
    % w=masse*u*v
    w=mtimeslike(SM.masse,u,@SEPtimesSEP_FH_MULTM,@SEPMATRIX)*v;
    % w=metric\w
    if SM.metric.m==1
        w=mtimeslike(SM.metric,w,@mldivide);
    elseif nargin>3 && isa(varargin,'SEPSOLVER')
        w=solve(SM.metric,varargin{1});
    else
        error('Donner un SEPSOLVER en argument')
    end
    % Demande de recompression de donnees
    if nargin>3, w=multisvd(w,varargin{:});end
else
    error('Non sens')
end



