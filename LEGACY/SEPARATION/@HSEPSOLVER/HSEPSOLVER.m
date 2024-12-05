function L = HSEPSOLVER(varargin)
% function L = HSEPSOLVER(dim,varargin)


% Les trucs qui servent � rien pour l'instant :
param = PARAMETERS(varargin{:});
% param = setdefaultparam(param,'type','alterne');
param = setdefaultparam(param,'adjoint',0);
% param = setdefaultparam(param,'adjointtype',0);
param = setdefaultparam(param,'metric',[]);
param = setdefaultparam(param,'inittype','rand');
% param = setdefaultparam(param,'restartifnotconverged',0);
% param = setdefaultparam(param,'errorindicator','none');
dim=5;% en attendant...
param = setdefaultparam(param,'updatedim',1:dim);
% param = setdefaultparam(param,'adjointdim',1:dim);
% param = setdefaultparam(param,'updateeps',1e-15);
% param = setdefaultparam(param,'ortho',0);
% param = setdefaultparam(param,'orthodim',1:dim);
% param = setdefaultparam(param,'updateadjoint',true);
% param = setdefaultparam(param,'deflationadjoint',true);

% param = setdefaultparam(param,'updatestep',1);
% param = setdefaultparam(param,'itercritupdate',1e-1);
% param = setdefaultparam(param,'reuse',false);
param = setdefaultparam(param,'righthandSD',false);
% param = setdefaultparam(param,'orthocrit',1e-8);
% param = setdefaultparam(param,'restart',0);
% param = setdefaultparam(param,'fullupdate',1);


% if isparamin(param,'reference') && ~isempty(getparam(param,'reference'))  
%   param = setparam(param,'errorindicator','reference');
% end

% param = setdefaultparam(param,'residual',false);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ceux qui se sont rajoutes en cours de route... %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param = setdefaultparam(param,'storeiter',0);
param = setdefaultparam(param,'updatestep',1);
param = setdefaultparam(param,'depth',0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ceux qui sont utiles (pour l'instant) dans HSM %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Le rapport aux arbres :
param = setdefaultparam(param,'node',1);
param = setdefaultparam(param,'tree',TREE());

% Si rien n'est precise (par defaut) :
nbnodes = length(getconnect(getparam(param,'tree')));
ONES    = ones(1,nbnodes);
ZEROS   = zeros(1,nbnodes);
DISP0   = ZEROS;
DISP0(1)=1;
% Les parametres importants :
param = setdefaultparam(param,'maxorder',  DISP0*12 + ONES* 3);
param = setdefaultparam(param,'maxiter',   ONES* 15);
param = setdefaultparam(param,'itercrit',  ONES* 1e-2);
param = setdefaultparam(param,'tol',       DISP0*1e-8 + (ONES-DISP0)* 5e-3);
% Les parametres moins importants :
param = setdefaultparam(param,'display' ,    DISP0);
param = setdefaultparam(param,'dyntol'  ,    ZEROS);
param = setdefaultparam(param,'dynitercrit', ZEROS);
param = setdefaultparam(param,'residual',    ZEROS);
ERRIND = cell(nbnodes,1);
ERRIND(:)={'none'};
param = setdefaultparam(param,'errorindicator',ERRIND);
% Les parametres d'update :
param = setdefaultparam(param,'alphaupdate', ZEROS);
param = setdefaultparam(param,'update',      ZEROS);
param = setdefaultparam(param,'updatetucker',ZEROS);

% Si juste la valeur est precisee :
% Les parametres importants :
if length(getparam(param,'maxorder'))    == 1 ,  param = setparam(param,'maxorder',   ONES*getparam(param,'maxorder')   );end
if length(getparam(param,'maxiter'))     == 1 ,  param = setparam(param,'maxiter',    ONES*getparam(param,'maxiter')    );end
if length(getparam(param,'itercrit'))    == 1 ,  param = setparam(param,'itercrit',   ONES*getparam(param,'itercrit')   );end
if length(getparam(param,'tol'))         == 1 ,  param = setparam(param,'tol',        ONES*getparam(param,'tol')        );end
% Les parametres moins importants :
if length(getparam(param,'display'))     == 1 ,  param = setparam(param,'display',    ONES*getparam(param,'display')    );end
if length(getparam(param,'dyntol'))      == 1 ,  param = setparam(param,'dyntol',     ONES*getparam(param,'dyntol')     );end
if length(getparam(param,'dynitercrit')) == 1 ,  param = setparam(param,'dynitercrit',ONES*getparam(param,'dynitercrit'));end
if length(getparam(param,'residual'))    == 1 ,  param = setparam(param,'residual',   ONES*getparam(param,'residual')   );end
if length(getparam(param,'errorindicator')) == 1
    ERRIND = cell(nbnodes,1);
    ERRIND(:)={getparam(param,'errorindicator')};
    param = setparam( param,'errorindicator',ERRIND );
end

param = setdefaultparam(param,'errorindicator',ERRIND);
% Les parametres d'update :
if length(getparam(param,'alphaupdate')) == 1 ,  param = setparam(param,'alphaupdate',ONES*getparam(param,'alphaupdate'));end
if length(getparam(param,'update'))      == 1 ,  param = setparam(param,'update',     ONES*getparam(param,'update')     );end
if length(getparam(param,'updatetucker'))== 1 ,  param = setparam(param,'updatetucker',ONES*getparam(param,'updatetucker'));end


% Ce qui peut (va) degager :
param = setdefaultparam(param,'startwith',1);
param = setdefaultparam(param,'start',0);

L=struct();
L = class(L,'HSEPSOLVER',param);



