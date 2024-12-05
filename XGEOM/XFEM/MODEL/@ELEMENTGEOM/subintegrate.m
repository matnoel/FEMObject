function Int = subintegrate(elem,xnode,subelem,subxnode,order,s,fun,varargin)

fun=fcnchk(fun);
if ~isa(subelem,'cell')
    subelem={subelem};
    subnode={subnode};
end
num = getnumber(elem);

%Int = zerosND([s,getnbelem(elem)]);
%for k=1:length(subelem)
%gauss=calc_gauss(subelem{k},order);

%numsub=getnumber(subelem{k});
%[temp,rep]=ismember(numsub,num);
%elemtemp = getelem(elem,rep);
%xnodetemp = xnode(:,:,rep);

%for l=1:gauss.nbgauss
% subdetJ=calc_detJ(subelem{k},subxnode{k},gauss.coord(l,:));
% xi = calc_x(subelem{k},subxnode{k},gauss.coord(l,:));
% x = calc_x(elemtemp,xnodetemp,xi);
% detJ=calc_detJ(elemtemp,xnodetemp,xi);
% f = fun(xi,elemtemp,xnodetemp,varargin{:}); 
% 
%Int(:,:,rep) = Int(:,:,rep) + gauss.w(l)*abs(subdetJ)*abs(detJ)*f;
%end
%end


Int = zerosND([s,getnbelem(elem)]);
for k=1:length(subelem)
gauss=calc_gauss(subelem{k},order);
numsub=getnumber(subelem{k});
[temp,rep]=ismember(numsub,num);
elemtemp = getelem(elem,rep);
xnodetemp = xnode(:,:,rep);
xgauss = gauss.coord;
wgauss = gauss.w;

subdetJ=calc_detJ(subelem{k},subxnode{k},xgauss);
xi = calc_x(subelem{k},subxnode{k},xgauss);
detJ=calc_detJ(elemtemp,xnodetemp,xi);
f = fun(xi,elemtemp,xnodetemp,varargin{:}); 

Int(:,:,rep) = Int(:,:,rep) + sum(wgauss*abs(subdetJ)*abs(detJ)*f,4);

end
