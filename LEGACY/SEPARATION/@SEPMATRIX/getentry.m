function a = getentry(v,i)

F = zeros(size(v.F));
  for k=1:v.m
            for nu=1:v.dim
    F(k,nu) = v.F{k,nu}(i(nu));
            end
  end
  
  a = prod(F,2)'*v.alpha(:);     
  