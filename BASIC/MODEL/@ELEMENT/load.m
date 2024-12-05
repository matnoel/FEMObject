function fe = load(elem,node,compo,fun,varargin)
% function fe = load(elem,node,compo,fun,varargin)
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
if isa(mat,'FOUR_ISOT')
    switch getdim(elem)
        case 1
            if isparam(mat,'S')
                S = evalparam(mat,'S',elem,xnode,gauss.coord); % cross-section area
                f = f*S;
            end
        case 2
            if isparam(mat,'DIM3')
                e = evalparam(mat,'DIM3',elem,xnode,gauss.coord); % thickness
                f = f*e;
            end
    end
end

if islocal(elem)
    P = calc_P(elem);
    N = N*P';
end

[repddl,repddlcompo] = findddl(elem.ddlnodedual,compo);

f = f(repddlcompo,:);

N = N(repddl,:);

fe = sum(gauss.w*abs(detJ)*N'*f,4);
