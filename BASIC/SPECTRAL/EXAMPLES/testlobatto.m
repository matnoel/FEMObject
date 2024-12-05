fe=1;

m = 3 ;
alpha = @(x) cos(pi*x*m);


if fe==0
    n = 15;
    La = POLYLAGRANGELOBATTO(n);
else
    nbelem = 10;
    n=10;
    La = POLYFELAGRANGE(linspace(-1,1,nbelem),n);
end
N = getnbpoints(La);
S = SPECTRALMODEL(1,La);
S = addbc(S,'dirichlet','all');
B = calc_matrix(S,1,1);
b = calc_vector(S,0,alpha);
%[A,B,C,b]=advec_conv_lobatto(La,alpha,'dirichlet');

uex = @(x) 1/pi^2/m^2*(cos(pi*x*m)-cos(pi*m));

u = B\b;
u=unfreevector(S,u);

figure(1)
clf
fplot(uex,[-1,1])
hold on
plot(u,S,'r')


%%
fe=1;
if fe==0
    n=20;
    La = POLYLAGRANGELOBATTO(n);
elseif fe==1
    nbelem=8;
    n=8;
    La = POLYFELAGRANGE(linspace(-1,1,nbelem),n);
end
S = SPECTRALMODEL(1,La);
casbc=2;
if casbc==1
    S = addbc(S,'dirichlet','right');
elseif casbc==2
    S = addbc(S,'periodic');
end
A = calc_matrix(S,0,0);
B = calc_matrix(S,1,1);
C = calc_matrix(S,0,1);
%[A,B,C,b]=advec_conv_lobatto(La,alpha,'periodic');

%
T = TIMEMODEL(0,.5,150);
k = .02;
v = -2;

N = DGTIMESOLVER(T,1);
u0 = @(x) exp(-x.^2*30);
%u0 = @(x) cos(pi*x);
%u0 = @(x) (1-x).*(1+x);
u0 = u0(getpoints(S));
u0 = u0(getddlfree(S));


ut = dsolve(N,zeros(size(A,1),length(N)),A,k*B+v*C,u0);
yax = [min(min(double(ut))),max(max(double(ut)))];

uref = ut;
if fe==0
    uref1=uref;
    La1 = La;
    S1=S;
else
    uref2=uref;
    La2=La;
    S2=S;
end

ut = setevolparam(ut,'ylim',yax);
ut = setevolparam(ut,'axison',true);

figure(1)

evol(ut,S)


%%

for i=1:length(uref1)
    u1 = uref1{i};
    u2 = uref2{i};
    figure(10)
    clf
    plot(u1,S1);
    hold on
    plot(u2,S2,'r')
    ylim(yax)
    pause(get(T,'t1')/length(N)/100)
end
uref = ut;

%%
x = getpoints(S);
ue=@(t) cos(pi*(x-v*t));
ue = @(t) cos(pi*x)*cos(pi*v*t) + sin(pi*x)*sin(pi*v*t)
tt=gettapprox(uref);

for i=1:length(uref)
    ti = tt(i);
    u = uref{i};
    figure(10)
    clf
    plot(u,S);
    hold on
    plot(ue(ti),S,'r')
    ylim(yax)
    pause(get(T,'t1')/length(N)/100)
end
%%
T = gettimemodel(N);

urad = TIMERADIALMATRIX(size(ut),T);

actu=1;

Md=12;
%Mt = getMmatrix(T);
%ft = zeros(size(A,1),length(T)) ;
%ft(:,1) = A*u0;
%ft = ft/Mt;
%ft = TIMEMATRIX(ft,T,size(u0));

ft = create_initial_condition(T,A*u0);


ftup=ft;

for i=1:Md
    l0 = rand(T);
    fprintf('--- fonction %d -----\n ',i)
    for j=1:5
        l=l0;
        AA = A *integratemtimes(l,diff(l)) + (k*B+v*C)*integratemtimes(l,l);
        bb = integratemtimes(ftup,l);
        V = AA\bb;
        V = V/norm(V);


        l = dsolve(N,V'*ftup,V'*A*V,V'*(k*B+v*C)*V);
        fprintf('iteration puissance : erreur %d\n',norm(l-l0)/norm(l));
        l0=l;
    end
    urad = urad+TIMERADIALMATRIX(V,size(ut),l);
    if actu==0
        ftup = ftup - (A*V)*diff(l) -((k*B+v*C)*V)*l;
    else

        L = getL(urad);

        itersub = 1;
        for kkk=1:max(itersub,1)
            AA = MULTIMATRIX(sparse(numel(A),numel(L)^2),size(A),[numel(L),numel(L)]);
            bb = MULTIMATRIX(sparse(size(A,1),numel(L)),[size(A,1),1],[numel(L),1]);
            for ii=1:numel(L)
                for jj=1:numel(L)
                    AA{ii,jj} = A *integratemtimes(L(ii),diff(L(jj))) + (k*B+v*C)*integratemtimes(L(ii),L(jj));
                    bb{ii} = integratemtimes(ft,L(ii));
                end
            end
            AA = assembleblock(AA);
            bb = assembleblock(bb);
            V = AA\bb;
            ftup = ft;
            V=MULTIMATRIX(V,size(ut),[1,numel(L)]);
            if itersub>0

                W =  double(V);
                L = dsolve(N,W'*ft,W'*A*W,W'*(k*B+v*C)*W);
            end
        end


        for ii=1:numel(L)
            ftup = ftup - (A*V{ii})*diff(L(ii)) -((k*B+v*C)*V{ii})*L(ii);
        end

        urad = TIMERADIALMATRIX(V,size(ut),L);


    end
%ftup = ftup - TIMEMATRIX((M*V)*(Dt*double(l)')' +(K*V)*(Mt*double(l)')',T,size(V));
    err(i)=norm(urad-uref)/norm(uref);

    fprintf('erreur : %d\n',err(i));
end



%%
for i=1:length(uref)
    u = getmatrixatstep(uref,i);
    ur = getmatrixatstep(urad,i);
    figure(10)
    clf
    plot(u,S);
    hold on
    plot(ur,S,'r')
    ylim(yax)
    pause(get(T,'t1')/length(N)/100)
end




%% STOCHASTIC
festo=  1;
cas=1;
if cas==1
    RV = RANDVARS(RVUNIFORM(-6,6));
elseif cas==2
    RV = RANDVARS(RVUNIFORM(0.01,0.1));
elseif cas==3
    RV = RANDVARS(RVUNIFORM(-1,1),RVUNIFORM(0.02,0.04));
elseif cas==4
    RV = RANDVARS(RVUNIFORM(-2,0),RVUNIFORM(0.02,0.04));
end
if festo
    X = PCMODEL(RV,'order',3*ones(1,getm(RV)),'fedim',1,'femesh',{10},'typebase',2);
else
    X = PCMODEL(RV,'order',5,'pcg');
end
PC = getPC(X);
if cas==1
    k=0.05;
    v = X{1};
elseif cas==2
    k=X{1};
    v = -2;
elseif cas==3
    k=X{2};
    v = X{1};
elseif cas==4
    k=X{2};
    v = X{1};
end
Mt = getMmatrix(N);
Dt=getDmatrix(N);
ft = zeros(size(A,1),length(N)) ;
ft(:,1) = A*u0;
ft = ft/Mt;
fsto=TIMEMATRIX(ft,N,size(u0));
Msto = A;
Ksto = B*k + C*v ;

%% calcul incremental de reference
N=setparam(N,'display',true);
[qpc,resultpcg] = dsolve(N,fsto,Msto,Ksto,[],@(A,b) cgs(A,b,1e-6));
%%
qpctemp = PCMATRIX(qpc);
qpcsd = spectral_decomposition(qpctemp,'tol',5e-3,'display',...
    'pfixmax',3,'pfixtol',1e-12,'nbfoncmax',30);
qpcsd = PCTIMEMATRIX(qpcsd,N);

%% Calcul spectral arnoldi global
G = GSDSOLVER('tol',1e-5,'nbfoncmax',35,...
    'reuse',true,...
    'display',true,...
    'errorindicator','residual','restart',7,'nbfoncmaxsimul',20,...
    'type','arnoldi','update',true,'orthocrit',1e-10,'direct',true);
N=setparam(N,'display',false);
[qrada,resulta] = dsolve(G,N,fsto.*one(PC),Msto*one(PC),Ksto);



%norm(qrada-qpc)/norm(qpc)
%%

%%
ut = random(qrada);
%ut = getV(qrada,7);
%ut = random(qpc);
yax = [min(min(double(ut))),max(max(double(ut)))];
ut = setevolparam(ut,'ylim',yax);
ut = setevolparam(ut,'axison',true);

figure(1)

evol(ut,S)
%%

[utemp,r] = random(qpc);
utempr = randomeval(qrada,r);
for i=1:length(uref)
    u = getmatrixatstep(utemp,i);
    ur = getmatrixatstep(utempr,i);
    figure(10)
    clf
    plot(u,S);
    hold on
    plot(ur,S,'r')
    ylim(yax)
    pause(get(T,'t1')/length(N)/100)
end
