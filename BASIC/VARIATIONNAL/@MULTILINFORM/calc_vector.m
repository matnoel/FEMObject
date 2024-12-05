function b = calc_vector(a,S,varargin)
% function b = calc_vector(a,S,varargin)

if isa(S,'SEPMODEL')
    b = calc_vector(S,a);
    return
end

if nargin==2
    varargin = cell(1,2);
    varargin(:) = {[]};
end

rep = find(~isemptyarg(varargin(1:getn(a))));
for i=1:length(rep)
    % try
    varargin{rep(i)} = unfreevector(S,varargin{rep(i)});
    % end
end

if ~isempty(a.k) && ~isempty(a.pk)
    try
        a.k = unfreevector(S,a.k);
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

if isa(getselgroup(a),'char')
    sel = 1:getnbgroupelem(S);
else
    sel = getselgroup(a);
end

be = cell(1,getnbgroupelem(S));
for p=sel
    be{p} = eval_elem(a,getgroupelem(S,p),getnode(S),varargin{:});
    trans = size(be{p},1)==1 && size(be{p},2)>1;
    if trans
        be{p} = be{p}';
    end
end

if israndom(be)
    b = assemble_pcvectorelem(S,be,'selgroup',sel);
else
    b = assemble_vectorelem(S,be,'selgroup',sel);
end

b = getfact(a)*b;

if isboundary(a)
    P = calc_P(Sini,S);
    b = P'*b;
end

if isfree(a)
    b = freevector(S,b);
end
