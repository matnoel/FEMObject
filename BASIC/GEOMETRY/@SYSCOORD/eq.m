function r = eq(s1,s2)
% function r = eq(s1,s2)

if (norm(s1.base-s2.base)<eps) && (s1.dim==s2.dim) && cellstrcmp(s1.axis,s2.axis)
    r = 1;
else
    r = 0;
end