function N = Ntet4(xi)
% function N = Ntet4(xi)

if nargin==0
    % N = inline('[1-xi(:,1)-xi(:,2)-xi(:,3),xi(:,1),xi(:,2),xi(:,3)]','xi');
    N = @(xi) [1-xi(:,1)-xi(:,2)-xi(:,3),xi(:,1),xi(:,2),xi(:,3)];
else
    % N=zeros(size(xi,1),4);
    % N(:,1)=1-xi(:,1)-xi(:,2)-xi(:,3);
    % N(:,2)=xi(:,1);
    % N(:,3)=xi(:,2);
    % N(:,4)=xi(:,3);
    N = [1-xi(:,1)-xi(:,2)-xi(:,3),xi(:,1),xi(:,2),xi(:,3)];
end
