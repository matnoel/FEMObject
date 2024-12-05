function u=switchmulti(u)

if isa(u.value,'double')
s=u.s;
sm=u.sm;
u.s=sm;
u.sm = s ;
u.value = u.value';
else
error('pas programme')
end

