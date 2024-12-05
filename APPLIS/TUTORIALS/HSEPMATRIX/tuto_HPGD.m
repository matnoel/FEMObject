
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% TEST multiSVD/PGD en HSEPMATRIX %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Charger les donnees :
% On charge les donnees d'un systeme lineaire :
%           Au=b

ABS=open('AUBsep.mat');
A=ABS.Af;
b=ABS.b;

PGD = SEPSOLVER(getdim(A),'maxorder',40,'tol',1e-5,'update',0);
u=solve(A,b,PGD);
PGDREF = SEPSOLVER(getdim(A),'maxorder',40,'tol',1e-5,'update',1);
uREF=solve(A,b,PGD);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% multiSVD %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tout commence avec un arbre :
%  1 2 3 4 5 6 7 8
c=[0 1 2 2 1 5 5 5];
T=TREE(c);

% On transforme u en HSEPMATRIX (qui servira de sol de reference):
HU = HSEPMATRIX(uREF,T);

% Quelques parametre de resolution :
%           1    2    3 4 5    6 7 8
maxorder = [10   3    0 0 3    0 0 0];
tol      = [1e-5 1e-2 0 0 1e-2 0 0 0];
alphaup  = [1    0    0 0 0    0 0 0];

% La compression en HmultiSVD :
HUsvd   = multisvd(HU,'maxorder',maxorder,'tol',tol);
HUsvdAU = multisvd(HU,'maxorder',maxorder,'tol',tol,'alphaupdate',alphaup);

% Tracer les convergences :
figure(1)
semilogy(fast_dist(HU     ,HU,HUsvd.m),'r')
hold on
semilogy(fast_dist(HUsvd,HU))
semilogy(fast_dist(HUsvdAU,HU,HUsvd.m),'y')
legend({ 'Reference','HmultiSVD','HmultiSVD avec AlphaUpdate' })
hold off


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% HPGD  %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arbre et parametre (les memes que precedemment):
%  1 2 3 4 5 6 7 8
c=[0 1 2 2 1 5 5 5];
T=TREE(c);
%           1    2    3 4 5    6 7 8
maxorder = [10   3    0 0 3    0 0 0];
tol      = [1e-5 1e-2 0 0 1e-2 0 0 0];
alphaup  = [1    0    0 0 0    0 0 0 ];
hsepsolver   = HSEPSOLVER('tree',T,'maxorder',maxorder,'tol',tol,'display',1);
hsepsolverAU = HSEPSOLVER('tree',T,'maxorder',maxorder,'tol',tol,'display',1,'alphaupdate',alphaup);

HA=HSEPMATRIX(A,T);
Hb=HSEPMATRIX(b,T);
HUpgd   = solve(HA,Hb,hsepsolver  );
HUpgdAU = solve(HA,Hb,hsepsolverAU);

figure(2)
semilogy( fast_dist(u,uREF,HUpgd.m) , 'r' )              %% REFERENCE
hold on
semilogy( fast_dist(HUpgd  , HSEPMATRIX(uREF,T) ) )      %% HPGD
semilogy( fast_dist(HUpgdAU, HSEPMATRIX(uREF,T) ) ,'y')  %% HPGD avec AU
legend({ 'Reference','HPGD','HPGD avec AlphaUpdate' })
hold off

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% HPGD profondeur 2  %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arbre et parametre (les memes que precedemment):
%  1 2 3 4 5 6 7 8 9 10
c=[0 1 2 3 3 2 6 1 8 8];
T=TREE(c);
%           1    2    3    4    5    6    7    8    9    10]
maxorder = [10   3    3    0    0    1    0    3    0    0];
tol      = [1e-5 1e-2 1e-2 0    0    1    0    1e-2 0    0];
alphaup  = [1    1    1    0    0    1    0    1    0    0];
disp     = [1    0    0    0    0    0    0    0    0    0];

hsepsolver   = HSEPSOLVER('tree',T,'maxorder',maxorder,'tol',tol,'display',disp);
hsepsolverAU = HSEPSOLVER('tree',T,'maxorder',maxorder,'tol',tol,'display',disp,'alphaupdate',alphaup);

HA=HSEPMATRIX(A,T);
Hb=HSEPMATRIX(b,T);
HUpgd   = solve(HA,Hb,hsepsolver  );
HUpgdAU = solve(HA,Hb,hsepsolverAU);

figure(3)
semilogy( fast_dist(u,uREF,u.m-2) , 'r' )                %% REFERENCE
hold on
semilogy( fast_dist(HUpgd  , HSEPMATRIX(uREF,T) ) )      %% HPGD
semilogy( fast_dist(HUpgdAU, HSEPMATRIX(uREF,T) ) ,'y')  %% HPGD avec AU
legend({ 'Reference','HPGD','HPGD avec AlphaUpdate' })
hold off



