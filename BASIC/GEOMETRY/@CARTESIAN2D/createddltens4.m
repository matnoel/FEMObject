function ddl = createddltens4(S,name)
% function ddl = createddltens4(S,name)
% name : nom du vecteur
% S : systeme de coordonnes

ddl = {[name '11'],[name '22'],[name '33'],...
    [name '12'],[name '13'],[name '23'],...
    [name '21'],[name '31'],[name '32']};
