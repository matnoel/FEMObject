function apc = project(a,pc,varargin)
% function apc = project(a,pc,varargin)

if isempty(getnumber(a))
    if getM(pc)==1 && isa(pc,'POLYCHAOS')
        % warning('la variable aleatoire n''a pas de numero : ca passe pour un chaos de dimension 1');
        a=setnumber(a,getnumber(pc,1));
    else
        error('donner un numero a la variable aleatoire correspondant a la dimension stochastique')
    end
end

[ok,rep]=ismember(a,pc);

if ~ok
    error('la variable aleatoire n''est pas contenue dans le chaos')
end


H = RANDPOLYS(pc);
H=H{rep};
p=getorder(pc);
p=p(rep);

apc = decomppc(a,H,'order',p,varargin{:});

apc = project(apc,pc);
