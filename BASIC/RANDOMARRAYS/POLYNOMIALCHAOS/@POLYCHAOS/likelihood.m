function L = likelihood(PC,a,Xs,varargin)
% L = likelihood(PC,a,Xs,varargin)
% PC : POLYCHAOS
% a : coefficients sur le polychaos
% appel de likelihood(X,Xs,varargin{:})
% avec X : PCMATRIX(a,s,PC)  (s : taille du vecteur aleatoire)

Xc = PCMATRIX(a,[size(a,1),1],PC);
L = likelihood(Xc,Xs,varargin);


