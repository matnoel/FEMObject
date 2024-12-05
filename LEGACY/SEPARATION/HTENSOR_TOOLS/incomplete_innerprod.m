function M = incomplete_innerprod(x,y,p)
% function M = incomplete_innerprod(x,y,p)
% Completely inspired from innerprod of htensor
M = cell(x.nr_nodes, 1);

x_is_leaf = x.is_leaf;
x_is_not_parent = setdiff(x.nr_nodes:-1:1,p);
x_is_not_parent = sort(x_is_not_parent,'descend');

% Start at leaves, move up the levels, dont treat parents
for ii=x_is_not_parent
    if(x_is_leaf(ii))
        % M_t = Ux_t' * Uy_t
        M{ii} = full(x.U{ii}'*y.U{ii});
    else
        ii_left  = x.children(ii, 1);
        ii_right = x.children(ii, 2);
        
        % M_t = Bx_t' * (M_tx kron M_ty) * By_t
        % (interpreting Bx_t, By_t to be in matricized form)
        B_ = ttm(y.B{ii}, { M{ii_left}, M{ii_right} }, [1 2]);
        
        M{ii} = ttt(x.B{ii}, B_, [1 2], [1 2], 3, 3);
        
        % If there is a singleton dimension, M{ii} needs to be
        % reshaped:
        M{ii} = reshape(M{ii}, size(x.B{ii}, 3), size(y.B{ii}, 3));
        
        % Save memory
        M{ii_left} = []; M{ii_right} = [];
    end
end

end


