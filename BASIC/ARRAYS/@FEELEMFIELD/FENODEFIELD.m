function un = FENODEFIELD(u)
% function un = FENODEFIELD(u)

if ~strcmp(u.storage,'node')
    error('il faut un stockage node du FEELEMFIELD')
end

v = u.value;
nbnode = size(v{1},3);
s = size(u);
un = zeros(nbnode,s(1),s(2));

repk = cell(1,length(v));
vk = cell(1,length(v));
for k=1:length(v)
    vk{k} = double(v{k});
    repk{k} = find(sum(sum(abs(vk{k}),1),2));
end

rep=[];
for k=1:length(v)
    repk{k} = setdiff(repk{k},rep);
    for i=1:s(1)
        for j=1:s(2)
            vij = vk{k}(i,j,:);
            vij = vij(:);
            un(repk{k},i,j) = vij(repk{k});
        end
    end
    rep = union(rep,repk{k});
end

un = FENODEFIELD(un(:,:));
