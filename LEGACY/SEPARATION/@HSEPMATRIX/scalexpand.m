function Hx = scalexpand(H)

% Pense dans le cas scalaire uniquement (cf norme)
% Hx=0; 
% alpha = H.alpha;
% for r=1:H.m
%     exr=alpha(r);
%     for d=1:H.dim
%         exr = exr * scalexpand(H.F{r,d});
%     end
%     Hx=Hx+exr;
% end



% exr=zeros(H.m,H.dim);
% for r=1:H.m
%     for d=1:H.dim
%     exr(r,d)=scalexpand(H.F{r,d});
%     end
% end
% Hx = H.alpha*prod(exr,2);

exr=zeros(H.m,H.dim);
for i=1:H.m*H.dim
    exr(i)=scalexpand(H.F{i});
end
Hx = H.alpha*prod(exr,2);