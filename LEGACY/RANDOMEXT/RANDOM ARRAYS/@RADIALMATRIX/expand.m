function v=expand(u)

v=(double(u.MULTIMATRIX)*u.D)*double(u.L);
v=MULTIMATRIX(v,u.s);

