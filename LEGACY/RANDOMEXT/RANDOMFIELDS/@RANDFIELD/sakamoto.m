function RFPC=sakamoto(RF,p,m,varargin)
% function RFPC=sakamoto(RF,p,m,S)
% RF : RANDFIELD
% S : MODEL de type MAILLAGE
% m : ordre de la decomposition de Karhunen-Loeve
% p : ordre du chaos
% RFPC : DISCRANDFIELD

disp('---- Decomposition Hermite + Karhunen-Loeve (Sakamoto,Ghanem 2002) ----')
S=getclassinvarargin('MODEL',varargin); 
S=MODEL('FIELD',{'U'},S);
D=double(calc_massegeom(S));
x = getcoord(S.node);

n= size(x,1);
h=POLYHERMITE();
xi=RANDVAR(h);
pc=POLYCHAOS(h,p);
gauss=calc_gausspoints(pc,max(getorder(pc))+1);
siz = [gauss.nbgauss,getP(pc)+1,n];
RFxi = margicdf(RF,x,cdf(xi,gauss.coord));
temp1 = repmat(pc(gauss.coord),[1,1,n]);
temp2 = repmat(RFxi,[1,1,getP(pc)+1]);
temp2=permute(temp2,[2,3,1]);

RFpc = gauss.w * (reshape(temp1,siz(1),prod(siz(2:3))).* ...
           reshape(temp2,siz(1),prod(siz(2:3))));
RFpc = reshape(RFpc,siz(2),siz(3))';

disp('Compute correlation matrix ...')
R=correl(RF,x,x);
s=std(RF,x);
s = s*s';
R = R.*s;

disp('COmpute correlation of gaussian field')
Rg=zeros(n,n);
% fun = inline('polyval([fliplr(u1(2:end).*u2(2:end)),-R12],x)','x','u1','u2','R12');
fun = @(x,u1,u2,R12) polyval([fliplr(u1(2:end).*u2(2:end)),-R12],x);
for i=1:n
for j=i:n
pourcentage((i)*n+j,n^2);
%r=roots([fliplr(upc(2:end).^2),-Ru(i,j)]);
%r=r(find(abs(imag(r))<eps));
%r = r(find(abs(r)==min(abs(r)))); 
%R(i,j)=sign(r)*min(abs(r),1);

Rg(i,j)=fzero(@(x) fun(x,RFpc(i,:),RFpc(j,:),R(i,j)),[-1.2,1.2]);
Rg(j,i)=Rg(i,j);

end
end
% Renormalisation de la correlation du champ gaussien
Rg=diag(sqrt(1./diag(Rg)))*Rg*diag(sqrt(1./diag(Rg)));


opt.disp=0;
disp('Solve eigenvalue problem ...')
[V,L]=eigs(Rg*D,m+1,'LM',opt);
lam=diag(L);
[lam,ia]=sort(lam,'descend');
error = sqrt(lam(end)/sum(lam));
disp(['decomposition KL du champ gaussien : ordre ' num2str(m)])
fprintf('   ->  error estimate = %3d\n',error);
lam = lam(1:m);
L = diag(lam);
V=V(:,ia(1:m));
V=V/diag(sqrt(diag(V'*D*V))); % D-normalisation des V
phi=diag(sqrt(sum((V*sqrt(L)).^2,2)))\(V*sqrt(L)); 

PC = POLYCHAOS(m,p,1);
P=getP(PC);
alm = getindices(PC);
RFPC=zeros(n,P+1);
prodphi=zeros(n,P+1);

prodphi=repmat(phi(:,1),1,P+1).^repmat(alm(:,1)',n,1);
for k=2:m
prodphi=prodphi.*(repmat(phi(:,k),1,P+1).^repmat(alm(:,k)',n,1));
end


t =(sqrt(factorial(alm(:,end)))./...
     sqrt(prod(factorial(alm(:,1:end-1)),2))');

RFPC = repmat(t,n,1).*RFpc(:,alm(:,end)+1).*prodphi;
RFPC=real(RFPC);


disp('Decomposition spectrale de la decomposition sur le chaos')
disp('Reduction du nombre de fonctions spatiales')

R = RFPC(:,2:end)*RFPC(:,2:end)';
[V,L]=eigs(R,m+1,'LM',opt);
lam=diag(L);
[lam,ia]=sort(lam,'descend');
error = sqrt(lam(end)/sum(lam));
fprintf('   ->  error estimate = %3d\n',error);
lam = lam(1:m);
L = diag(lam);
V=V(:,ia(1:m));
Vrad{1}=FENODEFIELD(RFPC(:,1));
lrad{1}=one(PC);
lpc = lrad{1}.value;
for i=1:m
 Vrad{i+1}=FENODEFIELD(V(:,i));
 lrad{i+1}=PCARRAY(PC,1/sqrt(lam(i))*V(:,i)'*(RFPC-repmat([RFPC(:,1)],1,P+1)));
 lpc = [lpc;lrad{i+1}.value];
end
lpc=PCARRAY(PC,lpc);
L=diag([1;sqrt(lam)]);
RFPC = PCRADIAL(Vrad,lpc,L);
%RFPC=PCARRAY(PC,RFPC);

RFPC = DISCRANDFIELD(RF,RFPC,S);