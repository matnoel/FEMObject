function u=switchmulti(u)

s=u.s;
sm=u.sm;
u.s=sm;
u.sm = s ;
u.value = u.value';


