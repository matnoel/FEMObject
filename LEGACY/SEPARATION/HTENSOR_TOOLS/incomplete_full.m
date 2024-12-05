function U = incomplete_full(x,p)
% function U = incomplete_full(x,p)

x_is_leaf = x.is_leaf;
x_is_not_parent = setdiff(x.nr_nodes:-1:1,p);
x_is_not_parent = sort(x_is_not_parent,'descend');

U = cell(1, x.nr_nodes);

% Loop through all nodes, starting with the leaf nodes.
for ii=x_is_not_parent

    if(x_is_leaf(ii))
% Matrix U already known
    U{ii} = x.U{ii};
else
% Find child nodes
    ii_left  = x.children(ii, 1);
    ii_right = x.children(ii, 2);

    BUU = ttm(x.B{ii}, {U{ii_left}, U{ii_right}}, [1 2]);
    U{ii} = matricize(BUU, [1 2], 3, false);

% Clear variables to save memory
    clear BUU;
    U{ii_left} = [];
    U{ii_right} = [];
end 

end


end


