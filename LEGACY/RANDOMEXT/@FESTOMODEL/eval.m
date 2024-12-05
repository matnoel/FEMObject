function a=eval(apc,x)
  
  a=zeros(size(x));
  for k=1:getM(apc)
  a(:,k) = eval(apc.xipc{k}',x);          
  end

  