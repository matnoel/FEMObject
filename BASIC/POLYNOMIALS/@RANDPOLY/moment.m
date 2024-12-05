function s=moment(polynome,e,varargin)

switch size(e,2)
case 2 % s = E(h_i * h_j)
if (e(1)==e(2))
s=1; % poly normalise
else
s=0;
end

case 3 % s = E(h_i * h_j * h_k)
i=e(1);j=e(2);k=e(3);

        rep0=find(e==0);

        
if ~isempty(rep0)
   repn0=setdiff([1:3],rep0(1));
   repn0=e(repn0);

   s=(repn0(1)==repn0(2));

else

if nargin<=2
intxn=calc_intxn(polynome,[0:sum(e)]);
else
intxn=varargin{1};
end

coeffi=polycoeff(polynome,i);
coeffj=polycoeff(polynome,j);
coeffk=polycoeff(polynome,k);
    
s=0;
for ll=0:i+j+k
coco = 0;
for ii=0:i
for jj=0:j
for kk=0:k
    if ((ii+jj+kk)==ll)
    coco=coco+coeffi(ii+1)*coeffj(jj+1)*coeffk(kk+1);
    end
end
end
end
s=s+coco*intxn(ll+1);
end

end   
  

    otherwise
fprintf('\n !!!! Coefficient Non defini !!!! \n')
stop
end

return


