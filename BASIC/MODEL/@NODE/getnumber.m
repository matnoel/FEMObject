function v = getnumber(u,listenode)
% function v = getnumber(u,listenode)

switch nargin
    case 1
        v = u.number;
    case 2
        v = u.number(listenode);
end
