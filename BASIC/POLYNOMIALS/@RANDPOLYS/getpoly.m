function h=getpoly(H,k)
% function h=getpoly(H,k)

if nargin==1
    h=H.h;
else
    h=H.h{k};
end

