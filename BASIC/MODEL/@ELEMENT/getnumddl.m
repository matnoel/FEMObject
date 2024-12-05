function n=getnumddl(elem,liste,choix)

if nargin==1
    n=elem.numddl;
    
else
    if nargin==3 & strcmp(choix,'global')
        liste = getpos(elem,liste);
    end
    n=elem.numddl(liste,:);
end