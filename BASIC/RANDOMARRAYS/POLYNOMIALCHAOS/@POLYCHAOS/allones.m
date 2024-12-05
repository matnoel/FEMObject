function onepc=allones(PC)
PC=getPC(PC);

onepc = ones(length(PC),1);
onepc = PCMATRIX(onepc,[1,1],PC);