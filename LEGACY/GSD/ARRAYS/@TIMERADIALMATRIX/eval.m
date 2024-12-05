function a=eval(trad,x)

a=double(trad.V) * trad.D* eval(trad.L,x);

  if size(a,2)==1
  a=reshape(a,size(trad));
  elseif all(size(trad)==1)  
  a=reshape(a,1,size(a,2));    
  else
  a=MULTIMATRIX(a,size(trad)); 
  end
