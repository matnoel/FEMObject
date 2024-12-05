function RF = RANDFIELD(varargin)
% function r = RANDFIELD(RFMARGINAL,RFCORREL)
% Creation d'un champ stochastique
% Le champ stochastique est exprim� comme une
% fonction non lineaire d'un champ gaussien centr� r�duit, 
% not� gamma(x,o) o� o est l'evenement. 
% RFMARGINAL est la loi marginale souhait�e pour le champ stochastique r(x,o)
% Si F_x est la fonction de distribution marginale de r en x,
% le champ stochastique s'exprime r(x,o)=F_x^-1(H(gamma(x,o))), o� H est la
% fonction de distribution gaussienne centr�e r�duite
% RFCORREL : fonction de correlation de gamma(x,o)
%
% function r = RANDFIELD(RFMARGINAL,RFCORREL,'notransform')
% Si option = 'selfcorrel', RFCORREL est alors la fonction de correlation de r
RF.marginal = varargin{1};
RF.correl = varargin{2};
RF.selfcorrel = ischarin('selfcorrel',varargin);

RF = class(RF,'RANDFIELD');
