function [rep,ia] = isin(p1,p2)

if p1.M<=p2.M  
 ia = zeros(1,p1.M) ;
  scan2=1:p2.M;
  for k=1:p1.M
    for k2=1:p2.M   
    if  any(k2==scan2) && polycmp(p1.h{k},p2.h{k2})
      ia(k)=k2;
      scan2(k2)=0;
      break
    end
    end
  end

  if any(~ia)
     ia=[];
     rep=0;
 else
 rep=1;
 end
else
rep=0;
ia=[];
end


