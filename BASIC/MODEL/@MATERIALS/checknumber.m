function checknumber(u)
% function checknumber(u)
% signale une erreur si deux materiaux ont le meme numero

number = getnumber(u);

if length(unique([number{:}]))~=length([number{:}])
    error('The materials must have different numbers.')
end
