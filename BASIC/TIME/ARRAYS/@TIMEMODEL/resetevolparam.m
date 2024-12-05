function T = resetevolparam(T)
% valeurs par defaut pour les options d'affichage

T.evolparam = PLOTOPTIONS();
T.evolparam = setparam(T.evolparam,'makepause',true,'pausetime',1/get(gettimemodel(T),'nt'),'plotiter',false,'plottime',true);
