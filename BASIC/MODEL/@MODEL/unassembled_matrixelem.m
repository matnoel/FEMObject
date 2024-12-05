function M = unassembled_matrixelem(S,me,varargin)
% s taille de la matrice finale
liste = getcharin('selgroup',varargin,1:S.nbgroupelem);
s=[S.nbddl,S.nbddl];
i=[];j=[];val=[];e=[];
for p=liste
    [ip,jp,ep,valp]=matrixelem(S.groupelem{p},double(me{p}),s);
    i=[i;ip];
    j=[j;jp];
    e=[e;ep];
    val=[val;valp];
end

% SPARSETENSOR VERSION
M=SPARSETENSOR([e i j],val);
