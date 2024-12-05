function [PCM,S]=PCMODEL(S,varargin)
% function PCM=PCMODEL(S,varargin)
% cree un PCMODEL a partir des variables aleatoires du MODEL
%    -> meme arguments que RANDVARS/PCMODEL
% Pour les dimensions stochastiques correspondant aux levelsets, on
% utilise des elements finis. Le maillage stochastique a été calculé par la
% fonction lsrandomsplit quand le modele a ete finalise. Ce maillage est
% stoke sous la forme de polynomes POLYFE dans S.H
%
% !!! si les levelsets et les materials dependent des memes variables aleatoires
% ne pas oublier de les numeroter avant de creer les levelsets ou les
% materials
%
% function [PCM,Spc]=PCMODEL(S,varargin)
% les parametres aleatoires des levelsets et materials de S sont remplaces
% par leur decomposition sur le chaos. les parametres deja decomposes sont
% projetes sur le nouveau chaos
% 
%
% -> See also RANDVARS/PCMODEL

RVls = RANDVARS(S.ls);
mat = MATERIALS(S);
RVmat = RANDVARS(mat);
RV = RANDVARS(RVls,RVmat);

n = getnumber(S.H);
if ~ischarin('manudim',varargin) & ~isclassin('RANDPOLYS',varargin)
PCM = PCMODEL(RV,'manudim',n,S.H,varargin{:});
else
error('manudim peut creer un conflit : on prefere rien faire')    
end

if nargout==2
S.ls = project(randomset(S.ls,PCM,RV),getPC(PCM));
mat = project(randomset(mat,PCM,RV),getPC(PCM));
S = actualisematerials(S); 
end



