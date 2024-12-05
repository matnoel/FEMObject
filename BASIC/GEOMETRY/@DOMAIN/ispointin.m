function [rep,P] = ispointin(D,P)
% function [rep,P] = ispointin(D,P)

c = getcoord(P);

P1 = double(getcoord(D.P1));
P2 = double(getcoord(D.P2));

switch D.dim
    case 1
        rep = find(c(:,1)>=P1(1)-eps & c(:,1)<=P2(1)+eps);
    case 2
        rep = find(c(:,1)>=P1(1)-eps & c(:,1)<=P2(1)+eps & ...
            c(:,2)>=P1(2)-eps & c(:,2)<=P2(2)+eps);
    case 3
        rep = find(c(:,1)>=P1(1)-eps & c(:,1)<=P2(1)+eps & ...
            c(:,2)>=P1(2)-eps & c(:,2)<=P2(2)+eps & ...
            c(:,3)>=P1(3)-eps & c(:,3)<=P2(3)+eps);
end

if nargout==2
    P = P(rep);
end
