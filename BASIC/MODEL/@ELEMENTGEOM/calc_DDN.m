function [DDN,detJ]=calc_DDN(elem,xnode,xgauss)

DDNlocal=getDDN(elem,xgauss);
  
% matrice des derivees des fonctions de forme par rapport aux variables
% locales
   								  % [N1,xi   N2,xi  ...
   								  %  N1,eta  N2,eta ... ]                          
[detJ,J,Ji] = calc_detJ(elem,xnode,xgauss);

% matrice des derivees des fonctions de forme par rapport aux variables
% globales
   								  % [N1,x   N2,x  ...
   								  %  N1,y  N2,y ... ]
n = getdim(elem);
if n==2
    Ji3 = [Ji(1,1)*Ji(1,1),Ji(1,2)*Ji(1,2),2*Ji(1,2)*Ji(1,1);...
        Ji(2,1)*Ji(2,1),Ji(2,2)*Ji(2,2),2*Ji(2,2)*Ji(2,1);...
        Ji(1,1)*Ji(2,1),Ji(1,2)*Ji(2,2),Ji(1,2)+Ji(2,1)]';
    DDN = Ji3*DDNlocal;
    
else
    error('pas programme')
    
end