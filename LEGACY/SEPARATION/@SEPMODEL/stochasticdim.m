function sd=stochasticdim(SM)
sd=[];
for d=1:SM.dim
    if strcmp(SM.F{d}.type,'RAND')
        sd=[sd d];
    end
end