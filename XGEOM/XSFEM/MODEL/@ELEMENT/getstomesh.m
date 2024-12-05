function n=getstomesh(elem,liste,choix)

if nargin==1
    n=elem.stomesh;
    
else
    if nargin==3 & strcmp(choix,'global')
        liste = getpos(elem,liste);
    end
    n=elem.stomesh{liste};
end