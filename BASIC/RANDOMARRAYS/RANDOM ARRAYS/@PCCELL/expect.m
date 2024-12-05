function Z = expect(K,a,b)
% function Z = expect(K,a,b)
% calcul de E(K*a)
% calcul de E(K*a*b)

s=size(K);
masse=getmasse(K);

switch nargin
    case 2
        PC = getPC(a);
        a=double(a);
        un = double(one(PC));
        for k=1:length(masse)
            massea{k} = a*masse{k}*un';
        end
        massea = getmatrix(massea);
        
    case 3
        PC = getPC(a);
        a=double(a);
        b = double(b);
        for k=1:length(masse)
            massea{k} = a*masse{k}*b';
        end
        massea = getmatrix(massea);
        
end

Z = getmatrix(K.value);
Z = Z*massea';
Z=reshape(Z,s);

