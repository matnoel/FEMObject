function s=moment(polynome,e,varargin)

switch size(e,2)
case 2 % s = E(h_i * h_j)
i=e(1);j=e(2);

coeffi=polycoeff(polynome,i);
coeffj=polycoeff(polynome,j);
pi = length(coeffi)-1;
pj = length(coeffj)-1;
if nargin<=2
intxn=calc_intxn(polynome,[0:pi+pj]);
else
intxn=varargin{1};
end

s=0;
for ll=0:pi+pj
coco = 0;
for ii=0:pi
for jj=0:pj
    if ((ii+jj)==ll)
    coco=coco+coeffi(ii+1)*coeffj(jj+1);
    end
end
end
s=s+coco*intxn(ll+1);
end
   

case 3 % s = E(h_i * h_j * h_k)
i=e(1);j=e(2);k=e(3);

coeffi=polycoeff(polynome,i);
coeffj=polycoeff(polynome,j);
coeffk=polycoeff(polynome,k);
pi = length(coeffi)-1;
pj = length(coeffj)-1;
pk = length(coeffk)-1;
if nargin<=2
intxn=calc_intxn(polynome,[0:pi+pj+pk]);
else
intxn=varargin{1};
end


s=0;
for ll=0:pi+pj+pk
coco = 0;
for ii=0:pi
for jj=0:pj
for kk=0:pk
    if ((ii+jj+kk)==ll)
    coco=coco+coeffi(ii+1)*coeffj(jj+1)*coeffk(kk+1);
    end
end
end
end
s=s+coco*intxn(ll+1);
end
   
  

    otherwise
fprintf('\n !!!! Coefficient Non defini !!!! \n')
stop
end

return


