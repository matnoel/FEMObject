function GSD = updatelocalstosolver(GSD)
% function GSD = updatelocalstosolver(GSD)

l0 = getparam(GSD,'localstosolver');
l1 = getparam(GSD,'localstosolveriter');
l2 = getparam(GSD,'localstosolverupdate');

if isempty(l0) && (isempty(l1) || isempty(l2))
   return
   %error('rentrer parametre localstosolver ou (localstosolveriter et localstosolverupdate) ') 
end

if isempty(l0) && ~isempty(l2)
    GSD = setparam(GSD,'localstosolver',l2);
elseif isempty(l0) && ~isempty(l1)
    GSD = setparam(GSD,'localstosolver',l1);
end
   
if isempty(l1) && ~isempty(l0)
    GSD = getparam(GSD,'localstosolveriter',l0);
end
if isempty(l2) && ~isempty(l0)
    GSD = getparam(GSD,'localstosolverupdate',l0);
end

