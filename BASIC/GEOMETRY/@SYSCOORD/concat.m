function r = concat(s1,s2)
% function r = concat(s1,s2)

if (s1.dim==s2.dim) && cellstrcmp(s1.axis,s2.axis)
    r = s1;
    r.MYDOUBLEND = concat(s1.MYDOUBLEND,s2.MYDOUBLEND);
end
