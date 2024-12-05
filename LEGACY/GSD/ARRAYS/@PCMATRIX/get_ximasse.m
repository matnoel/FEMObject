function ximasse=getximasse(upc)

ximasse = upc.ximasse;

if isempty(ximasse)
 error('ximasse pas calculee : utiliser calc_ximasse')      
end

