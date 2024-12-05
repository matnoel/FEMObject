function sto = FESTOMODEL(a,n,p,varargin)
% function sto = FESTOMODEL(a,n,p)
% a: RANDVARS (elles vont etre decomposees sur le une base EF)
% n : nombre de subdivisions de [0,1]^M
% n : tableau de cellules n{i} : contient le decoupage de [0,1]
% n : RANDPOLYS

a = RANDVARS(a);
disp('---- CREATION OF STOCHASTIC MODEL -----')    
M = getM(a);

disp('-> Unidimensional chaos decomposition of random variables')
[xipc,h] = decompfe(a,n,p);

disp('-> Creation of multidimensional chaos')
typebase = getcharin('typebase',varargin,2);
PC = POLYCHAOS(h,p,'typebase',typebase);

PC = calc_masseuni(PC);
PC = calc_masse(PC);

sto = struct();

for i=1:M
    sto.xipc{i}=prolonge(xipc{i},PC,i);
end
disp(' ')


sto=class(sto,'FESTOMODEL',PC);
superiorto('POLYCHAOS')
