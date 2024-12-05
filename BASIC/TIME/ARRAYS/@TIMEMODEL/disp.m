function disp(u)
% function disp(u)

a.t = u.t;
a.t0 = u.t0;
a.t1 = u.t1;
a.nt = u.nt;
a.dt = u.dt;
a.uniform = u.uniform;

disp(' MESH : ')
disp(a)

disp(' APPOXIMATION PARAMETERS : ')
disp(u.approxparam)




