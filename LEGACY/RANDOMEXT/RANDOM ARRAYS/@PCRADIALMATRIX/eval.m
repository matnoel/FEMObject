function a=eval(pcr,x)

a=double(pcr.V) * pcr.D* eval(pcr.L,x);

  if size(a,2)==1
  a=reshape(a,size(pcr));
  elseif all(size(pcr)==1)  
  a=reshape(a,1,size(a,2));    
  else
  a=MULTIMATRIX(a,size(pcr)); 
  end
