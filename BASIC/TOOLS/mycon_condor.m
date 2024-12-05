function [output,error]=mycon_condor(isGradNeeded,J,x,opt)

output = opt.mycon(isGradNeeded,J,x);

error=0;
