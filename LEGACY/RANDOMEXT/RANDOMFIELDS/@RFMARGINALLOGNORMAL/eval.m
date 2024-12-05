function rv = eval(ma,x)
% function rv = eval(ma,x)
% evaluation de la loi marginale aux points de coordonnees x

mu = evalparam(ma.mu,x);
sigma = evalparam(ma.sigma,x) ; 
x0 = evalparam(ma.x0,x) ; 

if ma.selfstat
  rv = RVLOGNORMAL(mu,sigma,x0,'stat')  ;   
else
  rv = RVLOGNORMAL(mu,sigma,x0)     ;
end


function val = evalparam(p,varargin)

    if isa(p,'inline') | isa(p,'function_handle')
    val=p(varargin{:});  
    else
    val=p;
    end