function fe = loadPG(elem,node,compo,fun,varargin)
% function fe = loadPG(elem,node,compo,fun,varargin)
% compo : nom des composantes de la force
% fun : voir comment la definir dans la fonction surfload

numelem = getnumber(elem);
xnode = getcoord(node,getconnec(elem)');
gauss = calc_gauss(elem,'mass');

x = calc_x(elem,xnode,gauss.coord);
[N,detJ] = calc_N(elem,xnode,gauss.coord);

if isa(fun,'FENODEFIELD')
    f = N*localize(elem,getvalue(fun));
else
    f = fun(x,varargin{:});
    f = MYDOUBLEND(f);
end

mat = getmaterial(elem);
if isparam(mat,'stabilize') && getparam(mat,'stabilize')>0
    b = evalparam(mat,'b',elem,xnode,gauss.coord) ;
    B = calc_B(elem,xnode,gauss.coord);

    bxi = norm(b);
    bxi = full(bxi);
    b = full(b);
    he = 2/sum(abs(b'*B/bxi));
    if getparam(mat,'stabilize')==1
        k = full(evalparam(mat,'k',elem,xnode,gauss.coord)) ;
        Pe = (bxi*he/2/k);
        tau = he/2./bxi.*(1./tanh(Pe)-1./Pe);
    elseif getparam(mat,'stabilize')==2
        tau = he/2/bxi;
    end

    N = N + tau*(b'*B);

end

fe = sum(gauss.w*abs(detJ)*N'*f,4);

