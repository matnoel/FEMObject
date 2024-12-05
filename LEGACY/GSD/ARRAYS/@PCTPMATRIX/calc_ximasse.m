function x = calc_ximasse(x)

PC = x.POLYCHAOSTP;
x.ximasse=cell(1,getnbgroups(PC));
for i=1:getnbgroups(PC)
if isranddim(x,i)
x.ximasse{i} = getmassegroup(PC,i);     
x.ximasse{i} = getmatrix(multimtimes(x.phi{i}',x.ximasse{i}));    
else
x.ximasse{i} = x.phi{i};  
end
end
