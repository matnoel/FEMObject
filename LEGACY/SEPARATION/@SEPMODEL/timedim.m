function sd=timedim(SM)
sd=[];
for d=1:SM.dim
    if strcmp(SM.F{d}.type,'TIME')
        sd=[sd d];
    end
end