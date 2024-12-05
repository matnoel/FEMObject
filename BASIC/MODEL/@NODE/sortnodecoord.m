function [u,oldnumber] = sortnodecoord(u,k)


if nargin==1 && getindim(u)>1
    error('rentrer une dimension pour le tri')
elseif nargin==1 && getindim(u)==1
    k=1;
end

co = getcoord(u);
[co,I]=sortrows(co,k);
u.POINT = u.POINT(I);
oldnumber = u.number(I);
u.number = [1:length(u.number)]';