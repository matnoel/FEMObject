function a = randomeval(rv,x,RV)
% function a = randomeval(rv,x,RV)
% RV : RANDVARS (M variables aleatoires)
% x : double n-by-M  ou  1-by-M cell de n-by-1 double (realisations de RV)
% a : n-by-1 double
% on cherche dans RV la variable de meme dimension stochastique que rv
% a sont alors les realisations correspondantes pour rv
% -> utilisation de [ok,rep] = ismember(rv,RV) et x = x(:,rep) ou x=x{rep}
% -> puis transfer(RV{rep},rv,x)

if nargin<3
    error('preciser les variables aleatoires correspondant aux realisations');
end

a = transfer(RV,rv,x);

