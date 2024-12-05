function errvec = compareQPTensors(x,y)

assert(x.space.dim(end)==y.space.dim(end),...
    'Ranks along micro dimension must be equal') ;
assert(strcmp(class(x.space),class(y.space)),'TSpaces must be of same type')
isOperator = isa(x.space,'TSpaceOperators') ;
errFun = @(a,b) norm(full(a-b)) ;
xs = x.space.spaces ;
ys = y.space.spaces ;
xc = double(x.core) ;
yc = double(y.core) ;
errvec = zeros(x.space.dim(end),3) ;
for n = 1:x.space.dim(end)
    xM = mesoFactor(xc,xs,n) ;
    yM = mesoFactor(yc,ys,n) ;
    errvec(n,1) = errFun(xM,yM) ;
    if isOperator
        errvec(n,2) = errFun(xs{end}{n},ys{end}{n}) ;
        xM = kron(xM,xs{end}{n}) ;
        yM = kron(yM,ys{end}{n}) ;
    else
        errvec(n,2) = errFun(xs{end}(:,n),ys{end}(:,n)) ;
        xM = kron(xM,xs{end}(:,n)) ;
        yM = kron(yM,ys{end}(:,n)) ;
    end
    errvec(n,3) = errFun(xM,yM) ;
end

end

function M = mesoFactor(xc,xs,n)

if numel(xs) == 2
    ixc = find(xc(:,n)) ;
    if iscell(xs{1}) % Operator
        M = xc(ixc(1),n)*xs{1}{ixc(1)} ;
        for i = ixc(2:end)
            M = M + xc(i,n)*xs{1}{i} ;
        end
    else % Vector
        M = xc(ixc(1),n)*xs{1}(:,ixc(1)) ;
        for i = ixc(2:end)
            M = M + xc(i,n)*xs{1}(:,i) ;
        end
    end
else % order 3
    [ixc,jxc] = find(xc(:,:,n)) ;
    if iscell(xs{1}) % Operator
        M = xc(ixc(1),jxc(1),n)*kron(xs{2}{jxc(1)},xs{1}{ixc(1)}) ;
        for k = 2:numel(ixc)
            M = M + xc(ixc(k),jxc(k),n)*kron(xs{2}{jxc(k)},xs{1}{ixc(k)}) ;
        end
    else % Vector
        M = xc(ixc(1),jxc(1),n)*kron(xs{2}(:,jxc(1)),xs{1}(:,ixc(1))) ;
        for k = 2:numel(ixc)
            M = M + xc(ixc(k),jxc(k),n)*kron(xs{2}(:,jxc(k)),xs{1}(:,ixc(k))) ;
        end
    end
end

end