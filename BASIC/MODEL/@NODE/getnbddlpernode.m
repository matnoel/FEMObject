function nb = getnbddlpernode(node,num,choix)
% function nb = getnbddlpernode(node,num,choix)

if nargin==1
    n = getnbddlpergroup(node);
    nb = n(node.repnodeingroupddl(:,1));
else
    if nargin==2
        choix = 'global';
    end
    if strcmp(choix,'global')
        pos = getpos(node,num);
    else
        pos = num;
    end
    nb = getnbddlpernode(node);
    nb = nb(pos);
    nb = reshape(nb,size(num));
end
