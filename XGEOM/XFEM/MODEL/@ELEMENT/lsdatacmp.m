function rep = lsdatacmp(elem1,elem2)

rep = strcmp(elem1.lstype,elem2.lstype);
rep = rep & all(elem1.lsnumber==elem2.lsnumber);
rep = rep & (elem1.lsenrich==elem2.lsenrich);
if ~isempty(elem1.lsnature) && ~isempty(elem2.lsnature)
rep = rep & strcmp(elem1.lsnature,elem2.lsnature);
end

if rep && elem1.lsenrich
   rep = rep & all(getparam(elem1,'lsenrichtype')==getparam(elem2,'lsenrichtype')); 
end
if rep && strcmp(elem1.lsnature,'crack') & (strcmp(elem1.lstype,'bicut') || strcmp(elem1.lstype,'tip')); 
   rep = rep & all(getparam(elem1,'tipnumber')==getparam(elem2,'tipnumber')); 
end
