function [DN,detJ] = calc_DN(elem,xnode,xgauss)
% function [DN,detJ] = calc_DN(elem,xnode,xgauss)
  
% matrice des derivees des fonctions de forme par rapport aux variables
% locales
   								  % [ N1,xi   N2,xi  ...
   								  %   N1,eta  N2,eta ... ]
                           
[detJ,J,Ji,xnode,DNlocal] = calc_detJ(elem,xnode,xgauss);

% matrice des derivees des fonctions de forme par rapport aux variables
% globales
   								  % [ N1,x  N2,x ...
   								  %   N1,y  N2,y ... ]

DN = Ji*DNlocal;
