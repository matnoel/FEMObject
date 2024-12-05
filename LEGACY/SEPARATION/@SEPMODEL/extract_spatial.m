function SM_spatial = extract_spatial(SM)
% SM_spatial est un SEPMODEL qui ne contient que la partie
% spatiale de SM
SM_spatial=extractdim(SM,spatialdim(SM));