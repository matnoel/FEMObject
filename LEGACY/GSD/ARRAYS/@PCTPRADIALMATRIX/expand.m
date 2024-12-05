function Z = expand(z)

Z = PCTPMATRIXSUM(z.POLYCHAOSTP);
for i=1:getm(z)
   Z = Z + getV(z,i)*getL(z,i); 
end