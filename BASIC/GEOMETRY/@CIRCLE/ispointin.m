function [rep,P] = ispointin(C,P)
% function [rep,P] = ispointin(C,P)

c = getcoord(P);
switch C.indim
    case 2
        rep = find(sqrt((c(:,1)-C.cx).^2+(c(:,2)-C.cy).^2)<=C.r+eps);
    case 3
        rep = find(sqrt((c(:,1)-C.cx).^2+(c(:,2)-C.cy).^2+(c(:,3)-C.cz).^2)<=C.r+eps);
end

if nargout==2
    P = P(rep);
end
