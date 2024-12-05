function h = POLYWAVELETS(n,p,varargin)
% function h = POLYWAVELETS(n,p)
% base multi-ondeletttes sur [0,1]
% polynomes par morceaux orthonorm�s sur [0,1] pour la mesure 1
% p : degre du polynome (ne sert pas a grand chose)
% n : resolution

if nargin==0
    h=struct();
    param=[];
    domain=[];
    h.basis=[];
    h.mother=[];

else      
    h=struct();
    param.n=n;
    for i=1:n
        param.I = linspace(0,1,2^i+1);
    end
    param.I = [param.I(1:end-1)',param.I(2:end)'];
    if nargin==1 || isempty(p)
        p=0;
    end
    param.p=p;
    domain = [0,1];
    p0 = POLYCHAOS(POLYFE(1,p),p);
    p1 = POLYCHAOS(POLYFE(2,p),p);
    p0 = decompmatrix(p1,[],[],@(x) polyval(p0,x(:))');
    p1 = decompmatrix(p1,[],[],@(x) polyval(p1,x(:))');
    for i=0:p
        p1(2*i+1) = p1(2*i+1) - p1(2*i+2);
    end
    p1 = p1(1:2:2*(p+1)); 
    for i=1:p+1
        for k=1:p+1    
            p1(i) = p1(i) - prodscal(p1(i),p0(k))*p0(k);       
        end
        for k=1:i-1
            p1(i) = p1(i) - prodscal(p1(i),p1(k))*p1(k);    
        end
        p1(i) = p1(i)/norm(p1(i));
    end
    h.basis = p0;
    h.mother = p1;
end

hp = RANDPOLY('wavelets',param,domain);
h = class(h,'POLYWAVELETS',hp);
superiorto('RANDPOLY');
