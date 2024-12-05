function setfemobjectoptions(option,val)
% function setfemobjectoptions(option,val)
% definition de parametres du code
% exemples : 
% 'tolerancepoint'
% 'tolerancelevelset'
% 'nbsamplingsrandomsplit'
% 'path'
% ...


global femobjectoptions

femobjectoptions = setparam(femobjectoptions,option,val);

