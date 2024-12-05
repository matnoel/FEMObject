function [rep,P] = ispointin(S,P)

rep=[];
for i=1:S.nbgroupelem
    [repP,dummy] = ispointin(S.groupelem{i},S.node,P);
    if ~isempty(repP)
        rep = union(rep,repP);
    end
end

if nargout==2
    P = P(rep);
end
