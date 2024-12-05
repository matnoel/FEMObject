function be = eval_elem(a,elem1,elem2,node1,node2,V1,V2,varargin)
% function be = eval_elem(a,elem1,elem2,node1,node2,V1,V2,varargin)

xnode = cell(1,2);
xnode{1} = node1(elem1);
xnode{2} = node2(elem2);
elem = {elem1,elem2};
node = {node1,node2};

if ischarin('quadorder',varargin)
    gauss = calc_gauss(a,elem{:},getcharin('quadorder',varargin));
else
    gauss = calc_gauss(a,elem{:});
end
xgauss = gauss.coord;
p = getp(a);
pk = getpk(a);


N = cell(1,2);
DN = cell(1,2);
B = cell(1,2);
N{1} = calc_N(elem{1},xnode{1},xgauss);
N{2} = calc_N(elem{2},xnode{2},xgauss);
[B{1},detJ,DN{1}] = calc_DNvect(elem{1},xnode{1},xgauss);
[B{2},detJ2,DN{2}] = calc_DNvect(elem{2},xnode{2},xgauss);
% detJ = calc_detJ(elem1,xnode1,xgauss);
V = {V1,V2};

repempty = find([isempty(V{1}),isempty(V{2})]);

for i=1:2
    if isempty(V{i})
        if p(i)==0
            V{i} = N{i};
        elseif p(i)==1
            V{i} = B{i};
            if ~isempty(a.dir)
                V{i} = V{i}(a.dir,:);
            else
                V{i} = sum(V{i}(1:getnbddlpernode(elem{i}),:),1);
            end
        end
    else
        if p(i)==0
            if isa(V{i},'cell') && isa(V{i}{1},'function_handle')
                V{i} = V{i}{1}(xglob);
            else
                V{i} = N*localize(elem{i},V{i});
            end
        elseif p(i)==1
            if isa(V{i},'cell') && isa(V{i}{1},'function_handle')
                V{i} = V{i}{1}(xglob);
                V{i} = MYDOUBLEND(V{i});
            else
                V{i} = B{i}*localize(elem{i},V{i});
                if ~isempty(a.dir)
                    V{i} = V{i}(a.dir,:);
                else
                    V{i} = sum(V{i}(1:getnbddlpernode(elem{i}),:),1);
                end
            end
        end
    end
end

be = V{1}'*V{2};
be = sum(gauss.w*abs(detJ)*be,4);


if isempty(repempty)
    be = double(sum(be,3));
end







