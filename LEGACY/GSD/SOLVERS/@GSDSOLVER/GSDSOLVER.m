function L = GSDSOLVER(varargin)
% function L = GSDSOLVER('type',type)
% Generalized Spectral Decoposition 
% type : 'power' (Power-type algorithm par defaut) 
%        'powerupdated' (Power-type algorithm with updating)
%        'arnoldi' (Arnoldi-type algorithm)
%        'powersubspace' subspace method
% function L = GSDSOLVER('nbfoncmax',nbfoncmax,'tol',tol)
% nbfoncmax : nombre maxi de fonctions radiales  (10 par defaut) 
% tol : tolerance pour la resolution  (10e-18 par defaut)
%
% function L = GSDSOLVER('pfixmax',pfixmax,'pfixtol',pfixtol)
% pfixtol : tolerance sur le point fixe (5e-2 par defaut)
% pfixmax : nombre maxi d'iterations de point fixe (3 par defaut)
%
% function L = GSDSOLVER('direct',direct,'toliter',toliter)
% direct : si true, resolution directe des equations stochastiques
%             false, resolution iterative (cgs)
% toliter : tolerance de l'algorithme iteratif pour la 
%             resolution des equations stochastiques (tol/10 par defaut)
%
% function L = GSDSOLVER('reuse',reuse)
% reutilisation d'une solution precedente (phase d'initialisation)
% reuse : true ou false
%
% function L = GSDSOLVER('display',display)
% display = true ou flase -> affiche la convergence de l'algorithme
%
% function L = GSDSOLVER('errorindicator',errorindicator)
% methode de calcul de l'indicateur d'erreur
% errorindicator = 'residual' -> calcul du residu
%                = 'rayleigh' -> base sur le quotient de Rayleigh ou
%                associe en non-symmetrique
param = PARAMETERS(varargin{:});
param = setdefaultparam(param,'type','power');
param = setdefaultparam(param,'inittype','random');
param = setdefaultparam(param,'nbfoncmax',10);
param = setdefaultparam(param,'nbfoncmaxsimul',getparam(param,'nbfoncmax'));
param = setdefaultparam(param,'orthocrit',1e-8);
param = setdefaultparam(param,'restart',0);
param = setdefaultparam(param,'pfixmax',3);
param = setdefaultparam(param,'orthoduringpfix',0);
param = setdefaultparam(param,'pfixtol',5e-2);
param = setdefaultparam(param,'tol',1e-10);
param = setdefaultparam(param,'tolini',getparam(param,'tol'));
param = setdefaultparam(param,'direct',false);
param = setdefaultparam(param,'toliter',min(getparam(param,'tol')/100,1e-2));
param = setdefaultparam(param,'tolupdate',min(getparam(param,'tol')/100,1e-2));
param = setdefaultparam(param,'display',false);
param = setdefaultparam(param,'errorindicator','residual');
param = setdefaultparam(param,'update',false);
param = setdefaultparam(param,'updatefactor',false);
param = setdefaultparam(param,'finalupdate',false);
param = setdefaultparam(param,'reuse',false);
param = setdefaultparam(param,'finalunique',false);
param = setdefaultparam(param,'finalSD',false);
param = setdefaultparam(param,'righthandSD',false);
param = setdefaultparam(param,'finalSDfacttol',0.1);
param = setdefaultparam(param,'finaluniquefacttol',0.1);
param = setdefaultparam(param,'orthogonalizeLresidual',true);
param = setdefaultparam(param,'subspaceiteration',0);
param = setdefaultparam(param,'localstosolver',[]);
param = setdefaultparam(param,'localstosolverupdate',[]);
param = setdefaultparam(param,'localstosolveriter',[]);
param = setdefaultparam(param,'adjoint',0);
param = setdefaultparam(param,'saveiter',false);

checkparamvalue(param,'type','power','arnoldi','powersubspace','power_separated');
checkparamvalue(param,'inittype','random','allone','one');
checkparamvalue(param,'display',true,false);
checkparamvalue(param,'errorindicator','residual','rayleigh','none','reference');

L=struct();
L = class(L,'GSDSOLVER',param);
