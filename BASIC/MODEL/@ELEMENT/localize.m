function qe = localize(elem,q,varargin)
% function qe = localize(elem,q,varargin)


if israndom(q)
    qe = localizepc(elem,q);
elseif isempty(q)
    qe = q;
else
    q = full(double(q));
    nbelem = getnbelem(elem);
    if ischarin('scalar',varargin)
        connec = getconnec(elem);
        qe = reshape(q(connec',:),[size(connec,2) size(connec,1) size(q,2)]);
    else
        nbddl = getnbddl(elem);
        numddlelem = getnumddl(elem);
        qe = reshape(q(numddlelem',:),[nbddl nbelem size(q,2)]);
    end
    qe = permute(qe,[1 3 2]);
    qe = MYDOUBLEND(qe);
    
end