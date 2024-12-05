function [ok,i1] = isin(pc1,pc2)
% function [ok,i1] = isin(pc1,pc2)
% ok=1 si pc1 est inclus dans pc2
% i1 donne le reperage des indices de pc1 dans pc2
% pc1.indices = pc2.indices(i1,:);


i1=[];

ok = isin(pc1.RANDPOLYS,pc2.RANDPOLYS);

if ok
ok = ok & all(pc1.p<=pc2.p);
if ok    
    [temp,i1]=ismember(pc1.indices,pc2.indices,'rows');
    if any(~temp)
        ok=0;
    end
end
end

