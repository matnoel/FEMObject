function S = createlsmodelfaces(S)

for i=1:getnbfacets(S)
temp = LSMODEL(getfacet(S,i),S.ls);    
S = setfacet(S,i,temp);
end

for i=1:getnbridges(S)
temp = LSMODEL(getridge(S,i),S.ls);    
S = setridge(S,i,temp);
end

