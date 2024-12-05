function [elems,reps] = lsrandomgetelem(elem,ls,choix,node)

switch choix
    case 'cut'
temp = ~(lsissurely(elem,ls,@lsisin,node) | lsissurely(elem,ls,@lsisout,node));
%temp = lsispossibly(elem,ls,@lsiscut,node);
reps = find(temp);
elems = getelem(elem,find(temp));
    case 'in'
temp = lsissurely(elem,ls,@lsisin,node);
reps = find(temp);
elems = getelem(elem,find(temp));
    case 'out'
temp = lsissurely(elem,ls,@lsisout,node);
reps = find(temp);
elems = getelem(elem,find(temp));
    otherwise
        error('mauvais choix')
end
