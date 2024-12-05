function L = SEPARATEDSOLVER(varargin)
% function L = SEPARATEDSOLVER()

param = PARAMETERS(varargin{:});
param = setdefaultparam(param,'onebyone',false);
param = setdefaultparam(param,'inittype','random');
param = setdefaultparam(param,'nbfoncmax',100);
param = setdefaultparam(param,'nbfoncmin',1);
param = setdefaultparam(param,'pfixmax',30);
param = setdefaultparam(param,'pfixstagn',10);
param = setdefaultparam(param,'tol',1e-7);
param = setdefaultparam(param,'pfixtol',getparam(param,'tol'));
param = setdefaultparam(param,'display',false);
param = setdefaultparam(param,'errorindicator','residual');
param = setdefaultparam(param,'update',true);
param = setdefaultparam(param,'cyclic',true);

checkparamvalue(param,'errorindicator','residual','none','reference');
checkparamvalue(param,'inittype','random','allone','one');

L=struct();
L = class(L,'SEPARATEDSOLVER',param);
