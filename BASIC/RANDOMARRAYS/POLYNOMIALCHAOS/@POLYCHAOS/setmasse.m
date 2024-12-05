function PC=setmasse(PC,masse)
% function PC=setmasse(PC,masse)

if ~isempty(masse) && (~isa(masse,'MULTIMATRIX') || length(masse)~=length(getPC(PC)))
    error('la masse stochastique n''a pas la bonne taille')
end

PC.masse = masse;

