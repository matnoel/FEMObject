%% Table de connectivite
% Les arbres sont des graphs non-orienté, acyclique et connexe. Ils
% representent la facon dont sont regroupee hierarchiquement les
% variables.
% 
% Un simple vecteur suffit a le definir. La ieme cordonnee renvoit au
% noeud auquel elle refere :

% noeud :  1 2 3 4 5 6
connect = [0 1 1 1 1 1];
figure(1)
treeplot(connect);
% On comprend que le noeud 1 est racine (0), et que les autres sont
% rattaches a 1. On peut ainsi definir un arbre plus complique :

% noeud :  1 2 3 4 5 6 7 8
connect = [0 1 1 2 2 3 3 3];
figure(2)
treeplot(connect);

%% Classe TREE
% Cette classe repond aux besoins specifiques de l'utilisation des arbres
% dans le code femobject.
% Les feuilles correspondent aux variables : elles sont automatiquement
% numerotee afin de pouvoir les identifier.
connect = [0 1 2 2 2 5 5 1 8 9 9 8 12 12 12];
T=TREE(connect);
figure(3)
plot(T);

%% Extraire une branche de l'arbre :
ST1 = subtree(T,2);
ST2 = subtree(T,8);
figure(4)
subplot(2,2,[1 2])
plot(T)
subplot(2,2,3)
plot(ST1)
subplot(2,2,4)
plot(ST2)

%% Assemblage d'arbres :
% L'ordre d'assemblage a de l'importance :
AT1 = TREE(ST1,ST2);
AT2 = TREE(ST2,ST1);
figure(5)
subplot(1,2,1)
plot(AT1)
subplot(1,2,2)
plot(AT2)

%% IMPORTANT : quelques conventions d'usage
% Lorsque l'arbre est construit a la main (avec une table de connectivite),
% il faut respecter quelques conditions :
%
% 1 : le premier noeud est TOUJOURS racine (connect = [0 ...])
%
% 2 : les feuilles sont toujours ratachees a un noeud qui ne contient que
%     des feuilles. Par exemple :
connect = [0 1 2 2 1];
INCORRECT=TREE(connect);
connect = [0 1 2 2 1 5];
CORRECT=TREE(connect);
figure(6)
subplot(1,2,1)
plot(INCORRECT);
subplot(1,2,2)
plot(CORRECT);
%     Sur l'arbre incorrect, un solve_alterne (par exemple) va appeler une
%     routine d'une SEPMATRIX sur le noeud 5 qui ne sera pas une SEPMATRIX 
%     (un simple double en l'occurence). L'arbre correct va creer ainsi une
%     SEPMATRIX de dimension 1 sur le neoud 5, ce qui resout le probleme.
%
% 3 : respecter un ordre bien precis pour l'appelation des noeuds. En
%     effet, l'ordre des noeuds influence leur nimerotation :
connect = [0 1 1 2 3 5 5 3 5 8 2 13 2 13 8];
INCORRECT=TREE(connect);
connect = [0 1 2 2 2 5 5 1 8 9 9 8 12 12 12];
CORRECT=TREE(connect);
figure(7)
subplot(1,2,1)
plot(INCORRECT);
subplot(1,2,2)
plot(CORRECT);
%     La regle de numerotation est la suivante : numeroter de bas en haut,
%     puis de gauche a droite. C'est un regle qui permet d'eviter toute
%     confusion dans les choix des dimensions.













