function v = expandprodtsepsep(T,S)
% function v = expandprodtsepsep(T,S)


v=expand(T*TSEPMATRIX(truncate(S,1)));
for i=2:S.m
    v=v+expand(T*TSEPMATRIX(truncate(S,i)));
end

