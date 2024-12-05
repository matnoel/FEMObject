function rv = RANDVAR(h)
    param = getparam(h);
    rv = RVDISCRETEUNIFORM(param.Q);
end
