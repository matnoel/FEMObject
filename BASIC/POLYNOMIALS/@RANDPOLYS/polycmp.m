function rep = polycmp(p1,p2)

rep= (p1.M == p2.M) ;
if rep
  for k=1:p1.M
   rep = rep & polycmp(p1.h{k},p2.h{k});
  if ~rep
   break
  end
  end     
end


