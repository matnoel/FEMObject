function M = unassembled_trimatrixelem(S,me,varargin)
% s taille de la matrice finale
% SPARSETENSOR par defaut (normal)

% ATTENTION : si plusieurs selgroup,
%      l'ordre global des elements
%      est celui des elements par
%      selgroup, puis dans l'ordre des
%      selgroup.
liste = getcharin('selgroup',varargin,1:S.nbgroupelem);
s=[S.nbddl,S.nbddl,S.nbddl];
i=[];j=[];k=[];val=[];e=[];
EP=0;
for p=liste
    [ip,jp,kp,ep,valp]=trimatrixelem(S.groupelem{p},double(me{p}),s);
    i =[i;ip];
    j =[j;jp];
    k =[k;kp];
    e =[e;ep+EP];
    EP=EP+sum(S.repelemingroupelem(:,1)==p);
    val=[val;(valp)];
end

M=SPARSETENSOR([e i j k],val);
M.bornes=[S.nbelem s];