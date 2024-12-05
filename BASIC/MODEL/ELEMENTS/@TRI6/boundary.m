function [segbord,nodebord,xnodeseg,repelembord] = boundary(elem,node)
% function [segbord,nodebord,xnodeseg,repelembord] = boundary(elem,node)

co = getconnec(elem);

seg = zeros(0,3);
seg = [seg;co(:,[1,2,4]);co(:,[2,3,5]);co(:,[3,1,6])];
nseg = size(seg,1);
segsort = sort(seg,2);
[segu,repu,rep2] = unique(segsort,'rows');
repnu = setdiff([1:nseg],repu);
segnu = segsort(repnu,:);
repnu = find(ismember(segsort,segnu,'rows'));

segbord = seg;
segbord(repnu,:) = [];

nodebord = getnode(node,(unique(segbord)));

numelem = getnumber(elem);
numelem = repmat(numelem(:),3,1);
numelem(repnu) = [];

xnodeseg = nodelocalcoord(elem);
xnodeseg = xnodeseg([1,2,4,2,3,5,3,1,6],:);
xnodeseg = reshape(xnodeseg,[3,3,2]);
xnodeseg = permute(xnodeseg,[1,3,4,2]);
xnodeseg = repmat(xnodeseg,[1,1,getnbelem(elem),1]);
xnodeseg = reshape(xnodeseg,[3,2,getnbelem(elem)*3]);
xnodeseg(:,:,repnu) = [];

segbord = SEG3(nodebord,numelem,segbord,...
    'material',getmaterial(elem),'option','BORD');

segbord = update_lsdata(segbord,elem);
