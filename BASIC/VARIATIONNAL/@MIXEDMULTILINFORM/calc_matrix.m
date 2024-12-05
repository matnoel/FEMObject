function b = calc_matrix(a,varargin)
% function b = calc_matrix(a,varargin)

S = varargin(1:getn(a));
varargin = varargin(getn(a)+1:end);

if nargin==1+getn(a)
    varargin = cell(1,getn(a));
    varargin(:) = {[]};
end

repmat = isemptyarg(varargin(1:getn(a)));
rep = find(~repmat);
repmat = find(repmat);

for i=1:length(rep)
    if isfloat(varargin{rep(i)}) && numel(varargin{rep(i)})>1
        varargin{rep(i)} = unfreevector(S{rep(i)},varargin{rep(i)});
    elseif isfloat(varargin{rep(i)})
        varargin{rep(i)} = varargin{rep(i)}*ones(getnbddl(S{rep(i)}),1);
    end
end


% if isboundary(a)
%     Sini = S;
%     S = create_boundary(Sini);
%     if ~isempty(getboundaryobject(a))
%         S = intersect(S,getboundaryobject(a));
%     end
%     S = createddlnode(S,getnode(Sini));
% end


if isa(getselgroup(a),'char')
    sel = 1:getnbgroupelem(S{1});
else
    sel = getselgroup(a);
end


be = cell(1,getnbgroupelem(S{1}));
for p=sel
    elem = cell(1,length(S));
    node = elem;
    for k=1:length(S)
        elem{k} = getgroupelem(S{k},p);
        node{k} = getnode(S{k});
    end
    be{p} = eval_elem(a,elem{:},node{:},varargin{:});
end
b = assemble_mixedmatrixelem(S{repmat(1)},S{repmat(2)},be,'selgroup',sel);
b = getfact(a)*b;
% if isboundary(a)
%     P = calc_P(Sini,S);
%     b = P'*b*P;
% end
if isfree(a)
    b = freevector(S{repmat(1)},b,1);
    b = freevector(S{repmat(2)},b,2);
end
end

