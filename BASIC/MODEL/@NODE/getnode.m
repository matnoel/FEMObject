function u = getnode(u,listenode,choix)
% function v = getnode(u,listenode,choix)
% si choix = 'global' (par defaut), listenode : numeros des noeuds
% si choix = 'local', listenode : position des noeuds dans la table
if nargin==1
    return
end
if nargin==2
    choix='global';
end
switch choix
    case 'global'
        listenode=listenode(:);
        pos=getpos(u,listenode);
    case 'local'
        pos = listenode;
    otherwise
        error('mauvais choix')
end


u.POINT  = u.POINT(pos(:));
u.number = u.number(pos(:));
u.nbnode = numel(u.number);

if ~isempty(u.lsenrichnode)
    rep = find(ismember(u.lsenrichnode,u.number));
    if ~isempty(rep)
        u.lsenrichnode = u.lsenrichnode(rep);
        u.repnature = u.repnature(rep);
        u.lsnumber = u.lsnumber(rep);
        u.lsenrichtype = u.lsenrichtype(rep);
    end
end



