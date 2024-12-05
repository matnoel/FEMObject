function checknumber(u)
% function checknumber(u)
% signale une erreur si deux variables ont le meme numero

number = getnumber(u);

if length(unique([number{:}]))~=length([number{:}])
    error('les variables doivent avoir des numeros differents')
end
