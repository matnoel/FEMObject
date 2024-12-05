function ke = stab(mat,elem,xnode,xgauss,varargin)
% function ke = stab(mat,elem,xnode,xgauss,varargin)

if getparam(mat,'stabilize')==0
    mat = setparam(mat,'stabilize',1);
end
k = full(evalparam(mat,'k',elem,xnode,xgauss));
B = calc_B(elem,xnode,xgauss);
b = getparam(mat,'b');

if isa(b,'cell')
    bx = b{1};
    by = b{2};
    
    bx = evalparam(setparam(mat,'b',bx),'b',elem,xnode,xgauss);
    by = evalparam(setparam(mat,'b',by),'b',elem,xnode,xgauss);
    
    b = [bx;by];
else
    
    if nargin>4 && isa(varargin{1},'POLYCHAOS')
        b = evalparampc(mat,'b',varargin{1},elem,xnode,xgauss);
    else
        b = evalparam(mat,'b',elem,xnode,xgauss);
    end
    
end

bxi = norm(b);

if getparam(mat,'stabilize')==3
    he = 2/sum(abs(b'*B/bxi));
    tau = he*bxi;
    ke = tau*B'*B;
elseif getparam(mat,'stabilize')==4
    he = calc_h(elem,xnode);
    tau = he*bxi;
    ke = tau*B'*B;
else
    b = full(b);
    he = 2/sum(abs(b'*B/bxi));
    
    if getparam(mat,'stabilize')==1
        Pe = (bxi*he/2/k);
        tau = he/2./bxi.*(1./tanh(Pe)-1./Pe);
    elseif getparam(mat,'stabilize')==2
        tau = he/2/bxi;
    end
    ke = tau*(B'*b)*(b'*B);
end
