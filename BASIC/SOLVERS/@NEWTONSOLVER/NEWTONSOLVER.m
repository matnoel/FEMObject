function N = NEWTONSOLVER(varargin)
% function N = NEWTONSOLVER('type',type)
% type = 'tangent' or 'full' : tangent Newton (default)
% type = 'modified' or 'pseudo' : modified Newton
% type = 'constant' or 'manual' : constant Newton
%
% function N = NEWTONSOLVER('maxiter',maxiter,'tol',tol)
% maxiter : maximum number of iterations (100 by default)
% tol : stopping criterion (1e-10 by default)
%
% function N = NEWTONSOLVER('tolreact',tolreact,'tolstagn',tolstagn)
% tolreact : reaction criterium for reactualisation of tangent matrix
%            (1e-1 by default)
% tolstagn : stagnation criterium for reactualisation of tangent matrix
%            (tol/10 by default)

param = PARAMETERS(varargin{:});
param = setdefaultparam(param,'type','full');
param = setdefaultparam(param,'increment',true);
param = setdefaultparam(param,'maxiter',100);
param = setdefaultparam(param,'tol',1e-10);
param = setdefaultparam(param,'tolreact',1e-1);
param = setdefaultparam(param,'tolstagn',getparam(param,'tol')/10);
param = setdefaultparam(param,'display',true);
param = setdefaultparam(param,'stopini',true);
checkparamvalue(param,'type','tangent','full','modified','pseudo','constant','manual');
checkparamvalue(param,'increment',true,false);

N = struct();
N = class(N,'NEWTONSOLVER',param);
