function [b,varargout] = init_resolution(T,b,varargin)
% function [b,varargout] = init_resolution(T,b,varargin)
% T : TIMEMODEL
% initialisation de la resolution d'une equation d''evolution
% transforme si necessaire le second membre en TIMEMATRIX
% ou si c'est une TIMEMATRIX, verifie la compatibilité avec le TIMEMODEL 

if nargout~=nargin-1
    error('rentrer le meme nombre d''arguments de sortie que d''entree')
end

if ~istime(b)
    if israndom(b)
    b = PCTIMEMATRIX(b,T,[size(b,1),1]);
    else
    b = TIMEMATRIX(b,T,[size(b,1),1]);    
    end
else
    b = transfer(b,T);
end

if length(gettimemodel(b))~=length(T)
    error('le second membre doit correspondre avec le TIMEMODEL')
end


varargout = varargin;