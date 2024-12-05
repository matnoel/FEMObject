function [facebord,nodebord] = boundary(elem,node)
% function [facebord,nodebord] = boundary(elem,node)

co = getconnec(elem);
repface = [1 2 3 5 6 7; 1 2 4 5 9 8; 1 3 4 7 10 8; 2 3 4 6 10 9];

face = zeros(0,6);
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

facebord = TRI6(nodebord,[1:size(facebord,1)],facebord,...
    'material',getmaterial(elem),'option','BORD');

facebord = update_lsdata(facebord,elem);
