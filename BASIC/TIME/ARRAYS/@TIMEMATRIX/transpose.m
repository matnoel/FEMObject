function u = transpose(u)
if length(u.s)>2
 warning('transposition sur les 2 premieres dimensions');
end

 s1=u.s(1);
 u.s(1)=u.s(2);
 u.s(2)=s1;


if isa(u.value,'cell')
for i=1:length(u.value)
    u.value{i} = tranpose(u.value{i});
end
else
if all(u.s>1)
    error('programmer que pour les vecteurs')
end
end