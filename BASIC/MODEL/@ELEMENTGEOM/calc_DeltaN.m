function DeltaN=calc_DeltaN(elem,xnode,xgauss)


DDN = calc_DDN(elem,xnode,xgauss);
   								  %  N1,y  N2,y ... ]
DeltaN = DDN(1:getdim(elem),:);

DeltaN = sum(DeltaN,1);
