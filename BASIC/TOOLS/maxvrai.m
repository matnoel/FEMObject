function a=maxvrai(p)
% fexp realisations de f
% p ordre du chaos sur lequel on veut l'identification de f
xi = RVNORMAL();
PC = POLYCHAOS(xi,p);
f = PCMATRIX([1,1,1/6],[1,1],PC);
fexp = random(f,1000);

options=optimset;
options = optimset(options,'Display','iter');
options = optimset(options,'TolFun',1e-2);
options = optimset(options,'TolX',1e-2);
options = optimset(options,'MaxIter',50);

m1 = 8;
m2 = 10;
da1 = [0.7,1.6];
da2 = [0,1/3];
[A1,A2] = meshgrid(linspace(da1(1),da1(2),m1),linspace(da2(1),da2(2),m2));
for i=1:size(A1,1)
    for j=1:size(A1,2)
        fprintf('[%d,%d]\n',A1(i,j),A2(i,j))    
        LA(i,j)=vrais([mean(fexp),A1(i,j),A2(i,j)],fexp,PC);
    end
end
figure(1)
surf(A1,A2,LA)
figure(2)
plot(A2(:,1),LA')
figure(3)
plot(A1(1,:),LA')
keyboard
%a0 = [mean(fexp),std(fexp),zeros(1,length(PC)-2)]; 
a0 = [0,1];
keyboard
a=fminsearch(@(a) vrais(a,fexp,PC) , a0,options);


function La=vrais(a,fexp,PC)

ns=5e3;
f = PCMATRIX(a,[1,1],PC);
fs = random(f,ns,1);
pfs = ksdensity(fs,fexp);
La=-sum(log(pfs));

return
