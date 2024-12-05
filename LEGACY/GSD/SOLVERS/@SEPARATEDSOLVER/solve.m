function [u,result] = solve(S,A,b,tol)
% function [u,result] = solve(S,A,b,tol)

if nargin==4
    S = setparam(S,'tol',tol);
end

if getparam(S,'onebyone')
    [u,result] = solve_onebyone(S,A,b);
else
    [u,result] = solve_block(S,A,b);
end

