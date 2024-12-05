function [P,rep1,rep2] = intersect(P1,P2)
% function [P,rep1,rep2] = intersect(P1,P2)
% P points communs aux points P1 et P2 (a la tolerance eps pres)
% P=P1(rep1) et P=P2(rep2)
%
% See also POINT/unique, POINT/setdiff, POINT/ispointin

rep1 = [];
rep2 = [];
P1 = POINT(P1);
P2 = POINT(P2);
if numel(P1)==0
    P = P1;
elseif numel(P2)==0
    P = P2 ;
else
    prec = getfemobjectoptions('tolerancepoint');
    
    coord1 = double(getcoord(P1));
    if ndims(coord1)>=4
        warning('on ne gere que les MYDOUBLEND a 3 dimensions')
    end
    
    coord1 = permute(coord1(:,:,:),[3 2 1]);
    coord1arr = funprec(coord1,prec);
    coord2 = double(getcoord(P2));
    if ndims(coord2)>=4
        warning('on ne gere que les MYDOUBLEND a 3 dimensions')
    end
    coord2 = permute(coord2(:,:,:),[3 2 1]);
    coord2arr = funprec(coord2,prec);
    
    [rep1,rep2]=ismember(coord1arr,coord2arr,'rows');
    rep1 = find(rep1);
    rep2 = rep2(rep1);
    P = POINT(permute(coord1(rep1,:),[3 2 1]));
end


function coord = funprec(coord,prec)

coord = round(coord/prec)*prec;

return
