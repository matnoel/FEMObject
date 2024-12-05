function N = Nseg2(xi)
% function N = Nseg2(xi)

if nargin==0
    % N = inline('[(1-xi(:))/2,(1+xi(:))/2]','xi');
    N = @(xi) [(1-xi(:))/2,(1+xi(:))/2];
else
    N = [(1-xi(:))/2,(1+xi(:))/2];   
end
