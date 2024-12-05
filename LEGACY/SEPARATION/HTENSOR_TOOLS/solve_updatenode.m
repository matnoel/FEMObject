function u = solve_updatenode(A,b,u,t)
% function u = solve_updatenode(A,b,u,t)

if u.is_leaf(t)
    u = solve_updateleafnode(A,b,u,t);
else
    u = solve_updatenonleafnode(A,b,u,t);
end

end
