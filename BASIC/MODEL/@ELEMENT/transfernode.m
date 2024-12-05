function se = transfernode(elem,node,se)
% function se = transfernode(elem,node,se)
% Pour determiner la valeur en un noeud, on fait la moyenne sur les elements du support du noeud 

connec = getconnec(elem);

se = double(se);
sse = size(se);
% se = reshape(se,[sse(1)*sse(2),sse(3)*sse(4)])';
se = reshape(se,sse(1)*sse(2),[])';
% sen = zeros(getnbnode(node),sse(1)*sse(2));
[temp,conneclocal] = ismember(connec,getnumber(node));
nodeconum = accumarray(conneclocal(:),1,[getnbnode(node),1]);

rep = conneclocal(:);
rep = repmat(rep,1,sse(1)*sse(2)) + ...
    repmat(getnbnode(node)*[0:sse(1)*sse(2)-1],length(rep),1);

sen = accumarray(rep(:),se(:),[getnbnode(node)*sse(1)*sse(2),1]);

sen = reshape(sen,[getnbnode(node),sse(1),sse(2)]);
sen = MYDOUBLEND(permute(sen,[2 3 1]));
repnz = find(nodeconum);
nodeconum = MYDOUBLEND(reshape(nodeconum,[1 1 length(nodeconum)]));
sen(:,:,repnz) = sen(:,:,repnz)/nodeconum(:,:,repnz);
se = sen;

% Pour determiner la valeur en un noeud, on choisit un element au hasard
% dans le support du noeud
% sse = size(se);
% se = reshape(se,[sse(1),sse(2) , sse(3)*sse(4)]);
% sen = zeros(sse(1),sse(2),getnbnode(node));
% [temp,conneclocal] = ismember(connec,getnumber(S.node));
% nodeconum = accumarray(conneclocal(:),1);
% [cou,b] = unique(conneclocal(:));
% sen(:,:,cou) = se(:,:,b);
% se = sen;