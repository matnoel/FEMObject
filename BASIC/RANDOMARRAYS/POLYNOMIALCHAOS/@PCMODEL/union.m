function PCM = PCMODEL(PCM1,PCM2,varargin)
% function PCM = PCMODEL(PCM1,PCM2,'typebase',typebase)
% union de deux PCMODEL

PC = union(PCM1.POLYCHAOS,PCM2.POLYCHAOS,varargin{:});

X1 = project(PCM1.X,PC);
X2 = project(PCM2.X,PC);

[temp,rep1]=ismember(PCM1,PC);
[temp,rep2]=ismember(PCM2,PC);

X = zeros(getM(PC),1,PC);
X(rep1) = X1(:);
X(rep2) = X2(:);

PCM = PCMODEL(X);
