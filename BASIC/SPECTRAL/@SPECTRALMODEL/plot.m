function plot(u,S,varargin)

u = unfreevector(S,u);
if S.dim==1
    npts = getcharin('npts',varargin,1000);
else
    npts = getcharin('npts',varargin,60);
end
varargin = delcharin('npts',varargin);


if S.dim==1
    domain = getdomain(S.L);
    x = linspace(domain(1),domain(2),npts);
    ux = eval_sol(S,u,x);
    
    plot(x,ux,varargin{:});
else
    domainx = getdomain(S.L);
    domainy = getdomain(S.L);
    
    x = linspace(domainx(1),domainx(2),npts);
    y = linspace(domainy(1),domainy(2),npts);
    [X,Y]=meshgrid(x,y);
    ux = eval_sol(S,u,[X(:),Y(:)]);
    ux = reshape(ux,size(X));
    surf(X,Y,ux,varargin{:});
end



