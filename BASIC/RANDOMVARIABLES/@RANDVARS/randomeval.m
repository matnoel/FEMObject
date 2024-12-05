function a = randomeval(rv,x,RV)
% function a = randomeval(rv,x,RV)
%
% a : n-by-M double ou M est le nombre de variable aleatoire rv
% x : n-by-M' double (ou 1-by-M' cell of n-by-1 double)
% x sont les realisations des M' variables RV
% les dimensions stochastiques de rv doivent etre contenues dans celles de RV
%      -> utilisation de transfer(RV,rv,x)
%
% See also RANDVARS/transfer, RANDVAR/transfer, RANDVAR/randomeval

if nargin<3
    error('preciser les variables aleatoires correspondant aux realisations');
end

a = transfer(RV,rv,x);

