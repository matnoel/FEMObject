function plot(h,liste)

if length(liste)>1
    hold on
end
for k=1:length(liste)
    v = zeros(1,length(h));
    v(liste(k)+1)=1;
    plot(PCMATRIX(v,[1,1],h));
end