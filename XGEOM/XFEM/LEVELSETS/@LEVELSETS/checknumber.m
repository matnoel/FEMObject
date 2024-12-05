function checknumber(u)
% function checknumber(u)
% signale une erreur si deux levelset ont le meme numero

number = getnumber(u);

if length(unique([number{:}]))~=length([number{:}])
    error('les levelset doivent avoir des numeros differents')
end
