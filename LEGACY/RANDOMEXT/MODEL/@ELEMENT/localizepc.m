function [qe,L] = localizepc(elem,q)

if ~israndom(q)
    qe = localize(elem,q);
else   
    
if isa(q,'PCMATRIX')
 L = getPC(q);
 q = double(q);
elseif isa(q,'PCRADIALMATRIX')
 L = getD(q)*getL(q);
 q = double(getV(q));
end

q=full(q);

nbddl = getnbddl(elem);
nbelem = getnbelem(elem);
numddlelem = getnumddl(elem);

qe=reshape(q(numddlelem',:),[nbddl 1 nbelem 1 size(q,2)]);
qe=PCMYDOUBLEND(MYDOUBLEND(qe),L,5);
end

