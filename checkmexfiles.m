mexfunctions={'multiplyF','redux'};

global mexcompiled;

if all(cellfun(@exist,mexfunctions)==3)
    mexcompiled=true;
else
    mexcompiled=false;
    warning('MexFilesCompilation',...
        'MEX FILES ARE NOT COMPILED, DO IT FOR SPEED ENHANCEMENT')
end
