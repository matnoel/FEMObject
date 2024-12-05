function [facebord,nodebord] = boundary(elem,node)
% function [facebord,nodebord] = boundary(elem,node)

co = getconnec(elem);
repface = [1 2 3; 1 2 4; 1 3 4; 2 3 4];

face = zeros(0,3);
for k=1:4
    face = [face;co(:,repface(k,:))];
end
nface = size(face,1);
facesort = sort(face,2);
[faceu,repu,rep2] = unique(facesort,'rows');
repnu = setdiff([1:nface],repu);
facenu = facesort(repnu,:);
repnu = find(ismember(facesort,facenu,'rows'));

facebord = face;
facebord(repnu,:) = [];

nodebord = getnode(node,(unique(facebord)));

facebord = TRI3(nodebord,[1:size(facebord,1)],facebord,...
    'material',getmaterial(elem),'option','BORD');

facebord = update_lsdata(facebord,elem);
