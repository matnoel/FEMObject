function sd=spatialdim(SM)
sd=[];
for d=1:SM.dim
    if strcmp(SM.F{d}.type,'SPACE')
        sd=[sd d];
    end
end