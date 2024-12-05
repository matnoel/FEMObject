function ls = uminus(ls)
% function ls = uminus(ls)

if iseval(ls)
    ls.value = -ls.value;
else
    ls.sign = -ls.sign;
end