function [elem,reps] = lsgetelem(elem,ls,choix,node)
% function [elem,reps] = lsgetelem(elem,ls,choix,node)

fun=fcnchk(['lsis' choix]);

reps=find(fun(elem,ls,node));
elem = getelem(elem,reps);

