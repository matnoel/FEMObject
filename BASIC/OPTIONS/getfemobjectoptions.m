function s = getfemobjectoptions(option)

global femobjectoptions
if isempty(femobjectoptions)
    initfemobjectoptions
end
s = getparam(femobjectoptions,option);

