function S=setoption(S,option,groupelem)
% function S=setoption(S,option,groupelem)

if nargin==2
  groupelem = 1:S.nbgroupelem ;  
end

for i=groupelem(:)'
    S.groupelem{i} = setoption(S.groupelem{i},option);
end