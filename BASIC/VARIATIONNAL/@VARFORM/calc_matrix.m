function b = calc_matrix(a,S,varargin)
% function b = calc_matrix(a,S,varargin)

rep = find(~isemptyarg(varargin(1:getn(a))));
for i=1:length(rep)
    varargin{rep(i)} = unfreevector(S,varargin{rep(i)});
end

try
    if ~isempty(getk(a)) && ~isempty(getpk(a))
        a = setk(unfreevector(S,getk(a)));
    end
end


if isboundary(a)
    Sini = S;
    S = create_boundary(Sini);
    if ~isempty(getboundaryobject(a))
        S = intersect(S,getboundaryobject(a));
    end
    S = createddlnode(S,getnode(Sini));
end

if isa(a.selgroup,'char')
    sel = 1:getnbgroupelem(S);
else
    sel = a.selgroup;
end
be = cell(1,getnbgroupelem(S));
for p=sel
    be{p} = eval_elem(a,getgroupelem(S,p),getnode(S),varargin{:});
end
b = assemble_matrixelem(S,be,'selgroup',sel);

b = getfact(a)*b;

if isboundary(a)
    P = calc_P(Sini,S);
    b = P'*b*P;
end

if isfree(a)
    b = freematrix(S,b);
end


