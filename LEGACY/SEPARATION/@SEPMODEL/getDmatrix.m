function Dt=getDmatrix(SM)
%     function Dt=getDmatrix(SM)
%
% Retourne l'operateur Dt de derivation temporel du
% SEPMODEL SM.
% Attention : une seule dimension temporelle est necessaire.

Dt=getmetric(SM);
locT=timedim(SM);
if ~all(size(locT)==1)
    error('Une dimension spatiale SVP');
end
dt=getDmatrix(SM.F{locT}.model);
Dt=fundim(Dt,@(a)dt,locT);