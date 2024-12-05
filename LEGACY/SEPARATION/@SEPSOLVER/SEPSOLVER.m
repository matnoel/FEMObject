function L = SEPSOLVER(dim,varargin)
% function L = SEPSOLVER(dim,varargin)

if nargin==0
    L = SEPSOLVER(0);
elseif isclassin('SEPSOLVER',varargin)
    L = getclassin('SEPSOLVER',varargin);
else
    
    param = PARAMETERS(varargin{:});
    param = setdefaultparam(param,'display',true);
    param = setdefaultparam(param,'type','alterne');
    param = setdefaultparam(param,'adjoint',0);
    param = setdefaultparam(param,'adjointtype',0);
    param = setdefaultparam(param,'metric',[]);
    
    param = setdefaultparam(param,'inittype','rand');
    
    param = setdefaultparam(param,'maxiter',10);
    param = setdefaultparam(param,'itercrit',5e-2);
    param = setdefaultparam(param,'restartifnotconverged',0);
    
    param = setdefaultparam(param,'tol',1e-8);
    param = setdefaultparam(param,'display',false);
    param = setdefaultparam(param,'errorindicator','none');
    
    param = setdefaultparam(param,'updatedim',1:dim);
    param = setdefaultparam(param,'updatedimtest',0);
    param = setdefaultparam(param,'updatedimapprox',[]);
    param = setdefaultparam(param,'adjointdim',1:dim);
    
    param = setdefaultparam(param,'update',0);
    param = setdefaultparam(param,'updateeps',0);
    param = setdefaultparam(param,'ortho',0);
    param = setdefaultparam(param,'orthodim',1:dim);
    % if length(getparam(param,'updatedim'))==1 && getparam(param,'adjoint')==0
    %     param = setparam(param,'update',min(getparam(param,'update'),1));
    % end
    
    param = setdefaultparam(param,'updateadjoint',false);
    param = setdefaultparam(param,'deflationadjoint',true);
    param = setdefaultparam(param,'alphaupdate',false);
    param = setdefaultparam(param,'updatestep',1);
    param = setdefaultparam(param,'itercritupdate',1e-1);
    
    param = setdefaultparam(param,'reuse',false);
    param = setdefaultparam(param,'righthandSD',false);
    param = setdefaultparam(param,'righthandSDstep',1);
    
    param = setdefaultparam(param,'orthocrit',1e-8);
    param = setdefaultparam(param,'restart',0);
    param = setdefaultparam(param,'fullupdate',1);
    param = setdefaultparam(param,'updatetucker',0);
    param = setdefaultparam(param,'seltuckerindices',0);
    
    
    param = setdefaultparam(param,'errorindicator','none');
    if isparamin(param,'reference') && ~isempty(getparam(param,'reference'))
        param = setparam(param,'errorindicator','reference');
    end
    
    param = setdefaultparam(param,'residual',false);
    param = setdefaultparam(param,'storeiter',0);
    
    % Solver non lineaire
    param = setdefaultparam(param,'bNlSolver',false);
    param = setdefaultparam(param,'cmdmLinOp',{});
    param = setdefaultparam(param,'modM1d',{});
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    param = setdefaultparam(param,'depth',0);
    param = setdefaultparam(param,'node',0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    param = setdefaultparam(param,'nodemaxiter',10);
    param = setdefaultparam(param,'nodetol',1e-4);
    param = setdefaultparam(param,'updateleaves',0);
    
    L=struct();
    L = class(L,'SEPSOLVER',param);
    
    
    
end
