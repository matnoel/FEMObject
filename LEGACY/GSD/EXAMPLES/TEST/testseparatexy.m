D=DOMAIN(2);
P=getvertices(D);
n=100;
[S,X,Y]=mesh(D,n,n);

S=setmaterial(S,FOUR_ISOT('k',1));
S=final(S);
S=addcl(S,P{4},'T',0);
A = calc_rigi(S);
%%
fun = @(x,y) exp(-100*(X-Y).^2);

uref=fun(X,Y);

figure(1)
clf
plot(uref(:),S,'surface')


[V,L]=eig(uref*uref');
ep=1e-3;
rep=find(sqrt(1-cumsum(flipud(diag(L)))./sum(diag(L)))<ep);
fprintf('critere = %f , nombre de modes = %d\n',ep,rep(1))

