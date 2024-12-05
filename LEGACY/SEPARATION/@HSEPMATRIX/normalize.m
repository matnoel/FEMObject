function [h,n] = normalize(h)
n=norm(h);
h.alpha = h.alpha/n;


