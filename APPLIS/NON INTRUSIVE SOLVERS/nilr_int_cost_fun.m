function [int] = nilr_int_cost_fun(V,Lam,cost_fun,gp,Hs)
% function [int] = nilr_int_cost_fun(V,Lam,cost_fun,gp,Hs)
%
% Compute the numerical integration of p \mapsto cost_fun(p,u(p))
% with u = V Lam'
%
% -gp contains the Gauss weights and points
%
% -Hs contains the evaluations of the basis of the stochastic space
% at the Gauss points


Uz = V*(Hs*Lam)';

int = zeros(gp.nbgauss,1);

coord = gp.coord;
for i = 1:gp.nbgauss
    int(i) = cost_fun(coord(i,:),Uz(:,i));
end

int = gp.w*int;

end