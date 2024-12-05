function [P,rep1,rep2] = setdiff(P1,P2)
% function [P,rep1,rep2] = setdiff(P1,P2);
% supprime de P1 les points communs aux points de P2 (a la tolerance eps pres)
% P=P1(rep1)
% P1(rep2) sont les points supprimes
%
% See also POINT/unique, POINT/intersect, POINT/ispointin

prec = getfemobjectoptions('tolerancepoint');
coord1 = permute(double(getcoord(P1)),[3 2 1]);
coord1arr = funprec(coord1,prec);
coord2 = permute(double(getcoord(P2)),[3 2 1]);
coord2arr = funprec(coord2,prec);
[rep1,rep2] = find(ismember(coord1arr,coord2arr,'rows'));
rep2 = find(rep2);
rep1 = setdiff(1:size(coord1,1),rep1);

P = POINT(permute(coord1(rep1,:),[3 2 1]));


function coord = funprec(coord,prec)

coord = round(coord/prec)*prec;

return