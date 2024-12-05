function ke = advdual(mat,elem,xnode,xgauss,varargin)
% function ke = advdual(mat,elem,xnode,xgauss,varargin)

B = calc_B(elem,xnode,xgauss);
N = calc_N(elem,xnode,xgauss);

b = getparam(mat,'b');

if isa(b,'cell')
    bx=b{1};
    by=b{2};
    
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

ke = (b'*B)'*N;
