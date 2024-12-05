function Nv = getNv(elem,xi)
% function Nv = getNv(elem,xi)

if nargin==2
    Nv = zeros([1,4,sizeND(xi)]);
    Nv(1,1,:) = 1/4.*(1-xi).^2.*(2+xi);
    Nv(1,2,:) = 1/8.*(1+xi).*(1-xi).^2;
    Nv(1,3,:) = 1/4.*(1+xi).^2.*(2-xi);
    Nv(1,4,:) = 1/8.*(1+xi).^2.*(xi-1);
else
    % Nv = inline('[(1-xi)^2*(2+xi)/4,1*(1+xi)*(1-xi)^2/8,(1+xi)^2*(2-xi)/4,1*(1+xi)^2*(xi-1)/8]','xi');
    Nv = @(xi) [(1-xi)^2*(2+xi)/4,1*(1+xi)*(1-xi)^2/8,(1+xi)^2*(2-xi)/4,1*(1+xi)^2*(xi-1)/8];
end

Nv = rescale(Nv,elem);
