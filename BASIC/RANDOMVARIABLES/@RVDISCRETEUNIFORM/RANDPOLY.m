function p = RANDPOLY(rv)
    param = getparam(rv);
    p=DIRACFUNCTIONS(param.Q);
end
