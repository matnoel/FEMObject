function GAUSS = subintegrategauss(elem,xnode,subelem,subxnode,order,s,fun,varargin)


if ~isa(subelem,'cell')
    subelem={subelem};
    subnode={subnode};
end
num = getnumber(elem);

for k=1:length(subelem)
gauss=calc_gauss(subelem{k},order);
numsub=getnumber(subelem{k});
[temp,rep]=ismember(numsub,num);
elemtemp = getelem(elem,rep);
subdetJ=calc_detJ(subelem{k},subxnode{k},gauss.coord);
xi = calc_x(subelem{k},subxnode{k},gauss.coord);
GAUSS{k}.coord = xi ;
GAUSS{k}.w = gauss.w*abs(subdetJ);
GAUSS{k}.nbgauss = gauss.nbgauss;
GAUSS{k}.repelem = rep;

end
