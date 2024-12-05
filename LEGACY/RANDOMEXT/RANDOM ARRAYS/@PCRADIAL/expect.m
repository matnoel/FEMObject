function Z = expect(K,a,b)
% function Z = expect(K,a,b)
% calcul de E(K*a)
% calcul de E(K*a*b)

switch nargin
    case 2
        
        Z = K.V{1}* expect(K.L{1},a);
        for k=2:K.m
            Z = Z + K.V{k} * expect(K.L{k},a);
        end
        
    case 3
        
        Z = K.V{1}* expect(a*K.Lmasse{1},b);
        for k=2:K.m
            Z = Z + K.V{k} * expect(a*K.Lmasse{k},b);
        end
        
end

