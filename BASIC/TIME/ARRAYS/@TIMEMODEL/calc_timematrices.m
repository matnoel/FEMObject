function T = calc_timematrices(T)

T.Mt = getMmatrix(T);
T.Dt = getDmatrix(T,'basic');
T.Diff = T.Mt\T.Dt;

