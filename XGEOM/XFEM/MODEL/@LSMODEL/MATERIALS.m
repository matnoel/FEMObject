function mat = MATERIALS(S)
% function mat = MATERIALS(S)

mat = MATERIALS(MATERIALS(S.MODEL),MATERIALS(S.ls));
