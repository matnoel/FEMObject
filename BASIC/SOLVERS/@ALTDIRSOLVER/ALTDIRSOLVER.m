function L = ALTDIRSOLVER(varargin)
% function L = ALTDIRSOLVER()
%
% function L = ALTDIRSOLVER('maxiter',maxiter,'tol',tol,'tolsyslin',tolsyslin)
% maxiter : maximum number of iterations (100 by default)
% tol : stopping criterion (1e-10 by default)
% tolsyslin : tolerance sur la resolution des systemes lineaires
%

param = PARAMETERS(varargin{:});
param = setdefaultparam(param,'maxiter',100);
param = setdefaultparam(param,'tol',1e-10);
param = setdefaultparam(param,'tolstagn',getparam(param,'tol')/10);
param = setdefaultparam(param,'display',true);
param = setdefaultparam(param,'stopini',true);

L=struct();
L = class(L,'ALTDIRSOLVER',param);
