function nbddl=getnbddl(SM)
% function nbddl=getnbddl(SM)
% nbddl(d) est le nombre de degre de liberte dans la dimension d

nbddl=cellfun(@(u) u.nbddl,SM.F);