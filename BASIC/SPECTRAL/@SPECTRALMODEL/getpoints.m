function x = getpoints(S)
% function x = getpoints(S)
if S.dim==1
    x=S.grid.X;
else
    x = [S.grid.X(:),S.grid.Y(:)]; 
end