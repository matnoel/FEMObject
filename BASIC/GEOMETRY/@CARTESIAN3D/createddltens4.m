function ddl = createddltens4(S,name)
% function ddl = createddltens4(S,name)
% name : nom du vecteur
% S : systeme de coordonnes

ddl = {[name '11'],[name '22'],[name '33'],[name '44'],[name '55'],[name '66'],...
    [name '12'],[name '13'],[name '14'],[name '15'],[name '16'],...
    [name '23'],[name '24'],[name '25'],[name '26'],...
    [name '34'],[name '35'],[name '36'],...
    [name '45'],[name '46'],...
    [name '56'],...
    [name '21'],[name '31'],[name '41'],[name '51'],[name '61'],...
    [name '32'],[name '42'],[name '52'],[name '62'],...
    [name '43'],[name '53'],[name '63'],...
    [name '54'],[name '64'],...
    [name '65']};
