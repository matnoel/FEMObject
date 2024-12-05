function plotmode(DRF,liste,varargin)

n=ceil(sqrt(length(liste)));
m=ceil(length(liste)/n);

for i=1:length(liste)
subplot(m,n,i)
plot(FENODEFIELD(getV(DRF.PCR,liste(i))),DRF.S,varargin{:});
title(['Mode ' num2str(liste(i))])
end