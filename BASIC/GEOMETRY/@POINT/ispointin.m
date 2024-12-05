function [rep,P] = ispointin(PO,P)
% function [rep,P] = ispointin(PO,P)
% cherche les points de P correspondant a PO
% si nargout == 1 : P(rep) sont les points communs
% si nargout == 2 : P correspond aux points communs
%
% See also POINT/unique, POINT/intersect, POINT/setdiff

[P,rep2,rep] = intersect(POINT(PO),POINT(P));
