function N = Ntri3(xi)
% function N = Ntri3(xi)

if nargin==0
    % N = inline('[1-xi(:,1)-xi(:,2),xi(:,1),xi(:,2)]','xi');
    N = @(xi) [1-xi(:,1)-xi(:,2),xi(:,1),xi(:,2)];
else
    N=zeros(size(xi,1),3);
    % N(:,1)=1-xi(:,1)-xi(:,2);
    % N(:,2)=xi(:,1);
    % N(:,3)=xi(:,2);
    N = [1-xi(:,1)-xi(:,2),xi(:,1),xi(:,2)];
end

