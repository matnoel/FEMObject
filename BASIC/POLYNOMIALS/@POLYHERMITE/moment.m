function s=moment(polynome,e,varargin)
switch size(e,2)
case 2 % s = E(h_i * h_j)
if (e(1)==e(2))
s=1; % poly normalise
else
s=0;
end

case 3 % s = E(h_i * h_j * h_k)

t=e(1)+e(2)+e(3);
if mod(t,2)==0  
p=t/2;
 if (p>=e(1) & p>=e(2) & p>=e(3))
 s=sqrt(factorial(e(1))*factorial(e(2))*factorial(e(3)))/...
 (factorial(p-e(1))*factorial(p-e(2))*factorial(p-e(3))); %poly normalise
 %s=(factorial(e(1))*factorial(e(2))*factorial(e(3)))/...
 %(factorial(p-e(1))*factorial(p-e(2))*factorial(p-e(3)));
 else
 s=0;
 end
else
s=0;
end

    otherwise
fprintf('\n !!!! Coefficient Non defini !!!! \n')
stop
end

return


