function xi = calc_xi(elem,xnode,x,varargin)
% function xi = calc_xi(elem,xnode,x,varargin)

display_ = ischarin('display',varargin);

if isa(xnode,'NODE')
    xnode = getcoord(xnode,getconnec(elem)');
end

if getindim(elem)~=getdim(elem)
    sys = getbase(getsyscoord(elem));
    sys = getsyscoordlocal(elem);
    xnode = xnode*sys;
    x = x*sys;
end

x1 = xnode(1,:);
x2 = xnode(2,:);
x3 = xnode(3,:);
x4 = xnode(4,:);

x = MYDOUBLEND(x);
xik = zerosND(size(x));

if display_
    fprintf('Fixed-point iteration for finding points in reference elements\n');
end
maxiter = 100;
tol = 1e-2;
for k=1:maxiter
    DN = calc_DN(elem,xnode,xik);
    xk = calc_x(elem,xnode,xik);
    % if norm(xk-x)<eps
    %    break
    % end
    fk = xk - x;
    alp0 = distance(POINT(x1),POINT(x3));
    alp = sqrt(sum(fk.*fk,2))./alp0;
    if display_
        fprintf('Iteration #%3.d : error = %.2d\n',k,max(double(alp),[],3));
    end
    if all(double(alp)<tol)
        break
    end
    
    df = DN*xnode;
    xik = xik - fk/df;
end

xi = xik;
