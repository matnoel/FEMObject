function L = GSDTSOLVER(varargin)
% function L = GSDTSOLVER('type',type)
% Generalized Spectral Decoposition 
% type : 'power' (Power-type algorithm par defaut) 
%        'powerupdated' (Power-type algorithm with updating)
%        'arnoldi' (Arnoldi-type algorithm)
%
% function L = GSDTSOLVER('nfoncmax',nbfoncmax,'tol',tol)
% nbfoncmax : nombre maxi de fonctions radiales  (10 par defaut) 
% tol : tolerance pour la resolution  (10e-18 par defaut)
%
% function L = GSDTSOLVER('pfixmax',pfixmax,'pfixtol',pfixtol)
% pfixtol : tolerance sur le point fixe (5e-2 par defaut)
% pfixmax : nombre maxi d'iterations de point fixe (3 par defaut)
%
% function L = GSDTSOLVER('direct',direct,'toliter',toliter)
% direct : si true, resolution directe des equations stochastiques
%             false, resolution iterative (cgs)
% toliter : tolerance de l'algorithme iteratif pour la 
%             resolution des equations stochastiques (tol/10 par defaut)
%
% function L = GSDTSOLVER('reuse',reuse)
% reutilisation d'une solution precedente (phase d'initialisation)
% reuse : true ou false
%
% function L = GSDTSOLVER('display',display)
% display = true ou flase -> affiche la convergence de l'algorithme
%
% function L = GSDTSOLVER('errorindicator',errorindicator)
% methode de calcul de l'indicateur d'erreur
% errorindicator = 'residual' -> calcul du residu
%                = 'rayleigh' -> base sur le quotient de Rayleigh ou
%                associe en non-symmetrique
param = PARAMETERS(varargin{:});
param = setdefaultparam(param,'type','power');
param = setdefaultparam(param,'inittype','one');
param = setdefaultparam(param,'nbfoncmax',10);
param = setdefaultparam(param,'nbfoncmaxsimul',10);
param = setdefaultparam(param,'orthocrit',1e-8);
param = setdefaultparam(param,'restart',0);
param = setdefaultparam(param,'pfixmax',3);
param = setdefaultparam(param,'pfixtol',5e-2);
param = setdefaultparam(param,'tol',1e-10);
param = setdefaultparam(param,'tolini',getparam(param,'tol'));
param = setdefaultparam(param,'direct',false);
param = setdefaultparam(param,'toliter',min(getparam(param,'tol')/100,1e-2));
param = setdefaultparam(param,'display',false);
param = setdefaultparam(param,'errorindicator','residual');
param = setdefaultparam(param,'update',false);
param = setdefaultparam(param,'finalupdate',false);
param = setdefaultparam(param,'reuse',false);
param = setdefaultparam(param,'finalunique',false);
param = setdefaultparam(param,'finalSD',false);
param = setdefaultparam(param,'righthandSD',false);
param = setdefaultparam(param,'finalSDfacttol',0.1);
param = setdefaultparam(param,'finaluniquefacttol',0.1);
param = setdefaultparam(param,'orthogonalizeLresidual',true);
param = setdefaultparam(param,'subspaceiteration',0);


checkparamvalue(param,'type','power','arnoldi','powersubspace');
checkparamvalue(param,'errorindicator','residual','rayleigh','none','reference');
checkparamvalue(param,'inittype','power','arnoldi','random','one');

L=struct();
L = class(L,'GSDTSOLVER',param);
