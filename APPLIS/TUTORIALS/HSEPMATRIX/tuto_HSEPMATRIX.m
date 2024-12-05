%% 
% Quelques vecteurs par dimensions (6)
u1= ones(5,1);
u2= ones(5,1);
u3= ones(5,1);
u4= ones(5,1);
u5= ones(5,1);
u6= ones(5,1);

% Un premier assemblage :
S1=SEPMATRIX({u1,u2;...
    u1,u2});         % SM de rang 2
S2=SEPMATRIX({u3,u4,u5,u6;
u3,u4,u5,u6;
u3,u4,u5,u6});   % SM de rang 3
H1=HSEPMATRIX({S1,S2;...
    S1,S2});        % HSM de rang 2

% H1 = sum^2(   sum^2( 12 )   .   sum^3( 3456 )  )
figure(1)
plot(H1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% On a acces a quelques infos :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total rank : rang de la representation SM de H1 :
TotalRank = totalrank(H1)
% Dimension de la HSM
Dim = H1.dim
% Arbre de construction de H1
T1  = H1.tree
% Norme :
NH1 = norm(H1)

%% On peut transformer H1 sous un autre format :

% Sepmatrix (SM)
SH1=SEPMATRIX(H1);

% En une autre HSM :
c2=[0 1 2 2 2 1 6 6 6];
T2=TREE(c2);
H2=HSEPMATRIX(H1,T2);
figure(2)
plot(H2)
RH2=H2.m

% La transformation precedante change la forme du premier niveau : H1 est
% alors developpee en SM, puis chaque rang est transforme en HSM selon T2.
% On compte alors 12 rang.

% -> La transformation est intelligente : elle conserve une representation
% la plus compacte. Ainsi, l'arbre 3 ne change pas le premier niveau de
% decomposition : la transformation ne s'applique qu'aux branches, et on
% conserve un rang 2 pour la premiere decomposition.
c3=[0 1 2 2 1 5 6 6 5 9 9];
T3=TREE(c3);
H3=HSEPMATRIX(H1,T3);
figure(3)
plot(H3)
RH3=H3.m

%%















