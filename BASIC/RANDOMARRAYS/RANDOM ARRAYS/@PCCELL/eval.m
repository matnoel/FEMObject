function a=eval(apc,x)
  
  PC=apc.POLYCHAOS;
  H=PC(x);
  a=H(1)*apc.value{1};
  for i=1:getP(PC)
  a=a+H(i+1)*apc.value{i+1};    
  end
  