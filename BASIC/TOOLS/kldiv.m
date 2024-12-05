function d = kldiv(x,p,q)
% divergence de kullback leibner

dx = x(2:end)-x(1:end-1);
p = (p(2:end)+p(1:end-1))/2;
q = (q(2:end)+q(1:end-1))/2;

d = sum((log(p(:))-log(q(:))).*p(:).*dx(:));
keyboard
