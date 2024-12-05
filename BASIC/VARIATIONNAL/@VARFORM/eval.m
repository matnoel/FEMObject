function ae = eval(a,S,varargin)


rep = find(~isemptyarg(varargin(1:getn(a))));
for i=1:length(rep)
varargin{rep(i)} = unfreevector(S,varargin{rep(i)});
end

try
if ~isempty(getk(a)) && ~isempty(getpk(a))
a = setk(unfreevector(S,getk(a)));   
end
end

if a.boundary
    Sini = S;
    S = create_boundary(Sini);
    if ~isempty(a.boundaryobject)
    S = intersect(S,a.boundaryobject);    
    end
    S = createddlnode(S,getnode(Sini));
end

if isa(a.selgroup,'char')
  sel = 1:getnbgroupelem(S);
else
  sel=a.selgroup;
end

ae = 0;
for p=sel
ae = ae + ...
    eval_elem(a,getgroupelem(S,p),getnode(S),varargin{:});
end

ae=getfact(a)*ae;



