function R=setximasse(R,DLmasse)
if ~isempty(DLmasse) && (~isa(DLmasse,'MULTIMATRIX') || length(DLmasse)~=getM(R))
    error('la masse stochastique n''a pas la bonne taille')
end

R.DLmasse=DLmasse;

