function X = partialprodscal(X,Y,dim)
% function Z = partialprodscal(X,Y,dim)

xd=X.dim;
xm=X.m;

xf=X.F;
yf=Y.F;

for dd=1:xd
    if dd~=dim
        for m=1:xm
            % Fastest of the 3 ways
            xf{m,dd} = sum(sum(xf{m,dd}.*yf{dd})); 
            % X.F{m,dd} = trace(X.F{m,dd}'*Y.F{dd});
            % X.F{m,dd} = sum(dot(X.F{m,dd},Y.F{dd}));
        end
    end
end

X.F=xf;

