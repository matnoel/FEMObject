function M = removeemptyfaces(M)
% function M = removeemptyfaces(M)

rep=[];
for i=1:length(M.facets)
    if getnbelem(M.facets{i})==0 || getnbnode(M.facets{i})==0
        rep=[rep,i];
    end
end
M.facets(rep)=[];

rep=[];
for i=1:length(M.ridges)
    if getnbelem(M.ridges{i})==0 || getnbnode(M.ridges{i})==0
        rep=[rep,i];
    end
end
M.ridges(rep)=[];


rep=[];
for i=1:length(M.peaks)
    if getnbelem(M.peaks{i})==0 || getnbnode(M.peaks{i})==0
        rep=[rep,i];
    end
end
M.peaks(rep)=[];


for i=1:length(M.facets)
    M.facets{i} = removeemptyfaces(M.facets{i});
end

for i=1:length(M.ridges)
    M.ridges{i} = removeemptyfaces(M.ridges{i});
end

for i=1:length(M.peaks)
    M.peaks{i} = removeemptyfaces(M.peaks{i});
end





