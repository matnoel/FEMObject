function v = mpower(u,m)
% function v = mpower(u,m)

if all(size2D(u)==1)
    v = double(u);
    v = MYDOUBLEND(reshape(v(:).^m,[1,1,sizeND(u)]));
elseif m==0
    v = repmat(eye(size2D(u)),[1,1,sizeND(u)]);
elseif isinteger(m) && m>0
    v = u;
    for i=2:m
        v = v*u;
    end
else
    error('mpower not defined for ND-array')
end