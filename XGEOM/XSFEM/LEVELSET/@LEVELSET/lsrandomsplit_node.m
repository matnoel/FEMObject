function [Ds,randh] = lsrandomsplit_node(ls,tol,varargin)
% function [Ds,randh] = lsrandomsplit_node(ls,tol,varargin)

D=getclassin('MODEL',varargin);


if ischarin('decomppc',varargin)
    PCM = PCMODEL(RANDVARS(ls),'order',4,'pcg');
    ls = project(ls,PCM);
    ls = lseval(ls,D);
end
X = RANDVARS(ls);



U = RANDVARS();
for i=1:getM(X)
    U{i}=RVUNIFORM(0,1);
end

D=getclassin('MODEL',varargin);


xitot = [0,1];
Iina = find_Iin(0,U,X,ls,D);
Iinb = find_Iin(1,U,X,ls,D);
xitot = recursivecut(xitot,0,1,Iina,Iinb,tol,U,X,ls,D);
xitot = sort(xitot);
randh = POLYFE(xitot);
Ds=D;

if ischarin('plot',varargin)
    plot(D)
    for i=1:length(xitot)
        x=transfer(U,X,xitot(i));
        lse = randomeval(ls,x,X);
        lse=lseval(lse,D);
        contourplot(lse);
    end
end

function xitot = recursivecut(xitot,a,b,Iina,Iinb,tol,varargin)
c=(a+b)/2;
Iinc = find_Iin(c,varargin{:});

if (b-a)<tol
    xitot = union(xitot,c);
else
    rep1 = comp_I(Iina,Iinc);
    rep2 = comp_I(Iinc,Iinb);
    
    if ~rep1
        xitot = recursivecut(xitot,a,c,Iina,Iinc,tol,varargin{:});
    end
    if ~rep2
        xitot = recursivecut(xitot,c,b,Iinc,Iinb,tol,varargin{:});
    end
    
end


return

function rep = comp_I(Ia,Ib)
rep = (length(Ia)==length(Ib)) && all(Ia==Ib) ;
return

function Iin = find_Iin(xi,U,X,ls,D)
x=transfer(U,X,xi);
ls = randomeval(ls,x,X);
ls=lseval(ls,D);
lss = sign(double(ls));
Iin = find(lss<0);
return
