function flag = check_reuse(GSD,ureuse)
% function flag = check_reuse(GSD,ureuse)
% verifie si la reutilisation peut ou doit etre faite

paramradial.reuse = getparam(GSD,'reuse');
if paramradial.reuse && ~isempty(ureuse)
   if isa(ureuse,'PCRADIALMATRIX')
        if length(ureuse)==0
            paramradial.reuse=false;
        end
    elseif ~isa(ureuse,'MULTIMATRIX')
        if isa(ureuse,'PCMATRIX')
           fprintf('  on ne peut reutiliser une PCMATRIX -> utiliser PCRADIALMATRIX ou MULTIMATRIX\n')
            paramradial.reuse=false;
        end
        if isa(ureuse,'double') & normest(ureuse)<eps
           fprintf('  vecteur de norme negligeable -> pas d''initialisation\n')
            paramradial.reuse=false; 
        end
   end
else
    paramradial.reuse = false;
end

flag= paramradial.reuse;