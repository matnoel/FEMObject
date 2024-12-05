function test = treecmp(T1,T2)
%   test = treecmp(T1,T2)
% Comparer le PREMIER etage des arbres
% test = 1 si les premieres decompositions sont les memes
% test = 0 si les premieres decompositions sont differentes

SN1 = find(T1.connect==1);
SN2 = find(T2.connect==1);
test=1;
if length(SN1) == length(SN2)
    for d=1:length(SN1)
        SV1=T1.Cvar{SN1(d)};
        SV2=T2.Cvar{SN2(d)};
        if length(SV1) == length(SV2)
        else
            test=0;
            return
        end
    end
else
    test=0;
    return
end



