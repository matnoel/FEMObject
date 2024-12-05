

%% VARIABLES ALEATOIRES

r = RVNORMAL(1,0.2);
pdfplot(r,'r')
random(r)

%% illustration d'un chaos generalise
PC = POLYCHAOS(RANDPOLYS(POLYHERMITE(),POLYLEGENDRE()),3)

getindices(PC)

%% decomposition d'une variable aleatoire sur le chaos
rl = RVLOGNORMAL(0.2,0.2,0,'stat');

rlpc = decomppc(rl,'order',3,'nbgauss',10);
getPC(rlpc) % pour obtenir le POLYCHAOS d'une PCMATRIX
figure(1000)
clf
pdfplot(rl)
hold on
pdfplot(rlpc,'r','nbs',1e5)

%% union de chaos

PC1 = getPC(rlpc);
PC1 = setnumber(PC1,1);
PC2 = POLYCHAOS(2,3);
PC2 = setnumber(PC2,2:3); % on associe des numeros aux dimensions stochastiques

PC = union(PC1,PC2);
PCnew = union(PC1,PC2);

rlpc2 = project(rlpc,PC);



%% CONSTRUCTION D'UN MODEL DE CHAOS
R = RANDVARS();
R{1} = RVNORMAL(1,0.2);
R{2} = RVUNIFORM(1,2);
X = PCMODEL(R,'order',5,'pcg');
PC = getPC(X);

pdfplot(X{1}+X{2},'nbs',1e5)


%% decomposition d'une fonction sur le chaos

fun = @(x) exp(3+.3*x);
PC = POLYCHAOS(POLYHERMITE(),1);
g = calc_gausspoints(PC,4)
fpc1 = decompfun(PC,13,[],fun);
PC = POLYCHAOS(POLYHERMITE(),3);
fpc2 = decompfun(PC,13,[],fun);
fr = fun(random(RVNORMAL(),10000,1));
figure(200)
clf
pdfsampleplot(fr,'bar','w');
hold on 
pdfplot(fpc1,'b','nbs',1e5)
pdfplot(fpc2,'g','nbs',1e5)

legend('reference tirages','p=1','p=3')




